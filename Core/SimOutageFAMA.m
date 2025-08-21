%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = SimOutageNakagami(Nsamples,gamma, Sigma, U, m) 
%
% Simulates the outage probability achieved by FAMA under Nakagami-m fading
% with correlation matrix Sigma. 
% 
% If U == 1 and gamma_avg ~= Inf, then it is a single user case and OP based 
% on SNR is simulated. Inteference is ignored.
%
% If U > 1 and gamma_avg == Inf, then noise is ignored and the OP based
% on SIR is simulated.
%
% If U > 1 and gamma_avg ~= Inf, then noise and interference are considered
% and OP based on SINR is simulated.
%
% Channel's complex Gaussian RVs are modeled acording with FAS "exact model" 
% presented in Khammassi's article [Eq. 10, 1], but with unitary variance of
% the i.i.d.normal RVs al and bl.
%
% Nakagami-m RVs are obtained from Gaussian RVs.
%
% For fast FAMA, it was used the inverse Nakagami-m CDF approximation [2].
%
% References:
% [1] M. Khammassi, A. Kammoun, and M.-S. Alouini, "A new analytical 
% approximation of the fluid antenna system channel,” IEEE Trans. 
% Wireless Commun., vol. 22, no. 12, pp. 8843–8858, Dec. 2023.
%
% [2] N. C. Beaulieu and C. Cheng, "Efficient Nakagami-m Fading Channel Simulation,"
% IEEE Trans. Veh. Techn., vol. 54, no. 2, pp. 413–424, Mar. 2005.
%  
% 
% -------------------- Parameters: ---------------------------------
%
% - Nsamples: number of Monte-Carlo simulations
% - gamma: SIR, SINR or SNR thresholds in linear scale (can be a vector).
% - gamma_avg: scalar, average received SNR (used for OP based on SNR and SINR)
% - U: number of users (scalar)
% - Sigma: correlation matrix of size NxN, where N is the number of ports.
%
% - m: severty of the fading in Nakagami-m distribution (integer)
%      for Fast and CFfast, m == {1,2,3,4,5,6,7,8,10} 
%
% - famatype: string that specifies the type of FAMA
%             If famatype == 'Fast', fast-FAMA
%             If famatype == 'Slow', slow-FAMA
%             If famatype == 'CFfast', cell-free fast-FAMA
%             If famatype == 'CFslow', cell-free slow-FAMA
%
%  - for free-cell MRT
%           d0 - ditance to target user
%           d - distance to interefence BS (all the same)
%           alphaCF - path loss exponent
%           Nant - number of antennas in the BS
%
% - pout: output vector (same size as gamma) with the outage probabilities.
%         pout(k) = P(SIR < gamma(k))
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = SimOutageFAMA(Nsamples,gamma, Sigma, U, m, famatype, gamma_avg, d0, d, alphaCF, Nant)



if strcmp(famatype, 'Fast') || strcmp(famatype, 'Slow') || strcmp(famatype, 'CFfast') || strcmp(famatype, 'CFslow')

    % Single user flag
    singleUser = (U == 1);

    % Pre-allocating
    pout = zeros(size(gamma));

    % For the sake of efficiency, only a few dominant eigenvectors and
    % eigenvalues can be computed, specially useful as N becomes large. If
    % desired, the full eigendecomposition can be computed always as
    % eig(Sigma).
    if numel(Sigma) > 1e5
        [V,Lambda] = eigs(Sigma,100);
    else
        [V,Lambda] = eig(Sigma);
    end

    % Sort eigenvalues in descending order
    [lambda, index] = sort(diag(Lambda),'descend');

    % Only the eigenvalues larger than a small tolerance are stored
    lambda = lambda(lambda>1e-5);

    % Re-arranges the corresponding eigenvectors
    V = V(:,index(1:length(lambda)));

    % Batch size to avoid RAM issues. Nsamples are generated in batches of
    % size batchsize
    batchsize = 1e3;

    % Generate variables by batches
    for batch = 1:batchsize:Nsamples

        % Single user
        if singleUser
            X = GenVariables_SingleUser(batchsize, V, lambda, m);
            if strcmp(famatype, 'Slow') || strcmp(famatype, 'Fast')
                sir = max(X ./ ((2*m)/(gamma_avg)),[],1);
            elseif strcmp(famatype, 'CFslow') || strcmp(famatype, 'CFfast')
                sir =  max( X ./ ((2*m*Nant)/(gamma_avg)),[],1);
            end
        % Multi-user
        else
            % Generate channels from desired user (X) and interference (Y)
            if strcmp(famatype, 'Slow')
                [Y, X] = GenVariables(batchsize, m*(U-1), V, lambda, m);
                Y = Y + ((2*m)/(gamma_avg));
            elseif strcmp(famatype, 'Fast')
                [Y, X] = GenVariablesFast(batchsize, U, V, lambda, m);
                Y = Y + ((2*m)/(gamma_avg));
            elseif strcmp(famatype, 'CFslow')
                [Y, X] = GenVariables(batchsize, m*(U-1), V, lambda, m*Nant);
                Y =  (d0/d)^alphaCF .* Y + ((2*m*Nant)/(gamma_avg));
            elseif strcmp(famatypeFast, 'CFfast')
                [Y, X] = GenVariablesFast(batchsize, U, V, lambda, m*Nant);
                Y =  (d0/d)^alphaCF .* Y + ((2*m*Nant)/(gamma_avg));
            end
           
            % Computes the SIR according to Eq. (24)
            sir = max(X./Y,[], 1);
        end
        % Computes outage probability
        for kg = 1:length(pout)
            pout(kg) = pout(kg) + sum(sir < gamma(kg));
        end
    end
    pout = pout/Nsamples;

else
    error('Type of FAMA must be indicaded: Fast, Slow, CFfast, CFslow')
end

%---------------------------------------------------------------------
%   Auxiliary function that generates correlated samples for the
%   channel from/to the desired user and the interference. The channel
%   is NAKAGAMI-m distributed with correlation matrix given by
%   Sigma = V'*diag(lambda)*V
%---------------------------------------------------------------------
    function [Y, X] = GenVariables(Nsamples, Util, V, lambda, m)

        N = size(V,1);
        %         disp ('size eigenvectors:')
        %         disp (N)

        % Pre-allocating
        % Squared Nakagami-m RVs
        X = zeros(N, Nsamples);
        Y = zeros(N, Nsamples);

        % Square-root of correlation matrix
        L = V*diag(sqrt(lambda));
        %         disp('size sqrt(correl matrix)')
        %         disp(size(L))

        % Number of independent Gaussians (rank of Sigma)
        Ngen = length(lambda);
        %         disp('length(lambda)')
        %         disp(Ngen)

        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Correlated samples of channel from/to desired user
            % Nakagami-m RV
            for l = 1:m
                aux = abs(L*(randn(Ngen, length(junk_index)) + ...
                    1i*randn(Ngen, length(junk_index)))/sqrt(2)).^2;
                X(:,junk_index) = aux + X(:,junk_index);
            end

            % Correlated samples of channel from/to interferers
            % and Nakagami-m RV
            for u = 1:(Util)
                aux = abs(L*(randn(Ngen, length(junk_index)) + ...
                    1i*randn(Ngen, length(junk_index)))/sqrt(2)).^2;
                Y(:,junk_index) = aux + Y(:,junk_index);
            end

        end

    end  % GenVariables


%---------------------------------------------------------------------
%   Auxiliary function that generates correlated samples for the
%   channel from/to the desired user and the interference. The channel
%   is NAKAGAMI-m distributed with correlation matrix given by
%   Sigma = V'*diag(lambda)*V
%   FAST FAMA
%---------------------------------------------------------------------

    function [Y, X] = GenVariablesFast(Nsamples, U, V, lambda, m)

        N = size(V,1);
        %         disp ('size eigenvectors:')
        %         disp (N)

        % Pre-allocating
        % Squared Nakagami-m RVs
        X = zeros(N, Nsamples);
        Y = zeros(N, Nsamples);
        % soma_aux = zeros(N, Nsamples);

        % Square-root of correlation matrix
        L = V*diag(sqrt(lambda));
        %         disp('size sqrt(correl matrix)')
        %         disp(size(L))

        % Number of independent Gaussians (rank of Sigma)
        Ngen = length(lambda);
        %         disp('length(lambda)')
        %         disp(Ngen)

        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Correlated samples of channel from/to desired user
            % Nakagami-m RV
            for l = 1:m
                aux = abs(L*(randn(Ngen, length(junk_index)) + ...
                    1i*randn(Ngen, length(junk_index)))/sqrt(2)).^2;
                X(:,junk_index) = aux + X(:,junk_index);
            end
            % Correlated samples of channel from interference users
            for u = 1:((U-1))
                gaux = L*( randn(Ngen, length(junk_index)) + 1i*randn(Ngen, length(junk_index)) )/sqrt(2);
                if m ~= 1
                    if m == 2
                        a1 = 0.1890; a2 = -0.0128; a3 = 0.2808; b1 = -0.0809; b2 = 0.0638;
                    elseif m == 3
                        a1 = 0.3472; a2 = -0.2145; a3 = 0.2626; b1 = -0.6779; b2 = 0.1690;
                    elseif m == 4
                        a1 = 0.4846; a2 = -0.4231; a3 = 0.2642; b1 = -0.9729; b2 = 0.2727;
                    elseif m == 5
                        a1 = 0.6023; a2 = -0.6238; a3 = 0.2789; b1 = -1.1798; b2 = 0.3732;
                    elseif m == 6
                        a1 = 0.7139; a2 = -0.8305; a3 = 0.3223; b1 = -1.3232; b2 = 0.4558;
                    elseif m == 7
                        a1 = 0.8167; a2 = -1.0244; a3 = 0.3761; b1 = -1.4233; b2 = 0.5192;
                    elseif m == 8
                        a1 = 0.9260; a2 = -1.2350; a3 = 0.4557; b1 = -1.4872; b2 = 0.5628;
                     elseif m == 10
                        a1 = 1.1088; a2 = -1.6095; a3 = 0.6015; b1 = -1.6046; b2 = 0.6488;                       
                    end
                    uniform = ( 1 - exp (- (abs(gaux).^2)./2) );
                    auxVar = (log(1./(1-uniform))).^(1/(2*m));
                    nkgNormalized = auxVar + ((a1.*auxVar + a2.*(auxVar.^2) + a3.*(auxVar.^3))./(1 + b1.*auxVar + b2.*(auxVar.^2)));
                    nkg = sqrt(2) * nkgNormalized ;
                    nkgComplex = nkg .* exp(1i.*angle(gaux));
                end
                if m == 1
                    nkgComplex = gaux;
                end
                Y(:,junk_index) = nkgComplex + Y(:,junk_index);
            end

            Y(:,junk_index) = abs(Y(:,junk_index)).^2;

        end

    end % function GenVariablesFast


%---------------------------------------------------------------------
%   Auxiliary function that generates correlated samples for the
%   channel from/to the desired user. The channel
%   is NAKAGAMI-m distributed with correlation matrix given by
%   Sigma = V'*diag(lambda)*V
%---------------------------------------------------------------------
    function [X] = GenVariables_SingleUser(Nsamples, V, lambda, m)

        N = size(V,1);

        X = zeros(N, Nsamples);

        L = V*diag(sqrt(lambda));

        Ngen = length(lambda);

        % Loop over junks to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Nakagami-m RV
            for l = 1:m
                aux = abs(L*(randn(Ngen, length(junk_index)) + ...
                    1i*randn(Ngen, length(junk_index)))).^2;
                X(:,junk_index) = aux + X(:,junk_index);
            end

        end

    end % function GenVariables_SingleUser

end

