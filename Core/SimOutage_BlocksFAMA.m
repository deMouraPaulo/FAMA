%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = SimOutage_BlocksNakagami(Nsamples,gamma, U, mu, L, m)
%
% Simulates the outage probability achieved by FAMA under Nakagami-m fading
% with block-diagonal correlation matrix. 
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
% 
%  ------------------ Parameters: ----------------------------------------
%
% - Nsamples: number of Monte-Carlo simulations;
% - gamma: SIR threshold (can be a vector);
% - gamma_avg: scalar, average received SNR (used for OP based on SNR and SINR);
% - U: number of users (scalar);
% - mu: correlation factor (scalar) within each block
%       obs.: mu = sqrt(delta), where delta is the power correlation coefient;
% - L: vector containing the size of each block;
% - m: Nakagami-m fading severity;
% - famatype: string that specifies the type of FAMA
%             If famatype == 'Fast', fast-FAMA;
%             If famatype == 'Slow', slow-FAMA;
%             If famatype == 'CFfast', cell-free fast-FAMA;
%             If famatype == 'CFslow', cell-free slow-FAMA;
%   
%   - flagmf: flag for approx. method of simulation for fast-FAMA
%       if = flagmf = 1 (true), approximated method (mf = 1, approx. method, Gaussian RVs)
%       if = flagmf = 0 (false), Nakagami-m RVs
%
%
%  - for free-cell FAMA with MRT precoding:
%       d0 - distance from serving BS to the target user;
%       d - distance from interfering BSs to the target use; 
%               (d is the same for all interfering BS);
%           alphaCF - path loss exponent;
%           Nant - number of antennas in the BS.
%
% - pout: output vector (same size as gamma) with the outage probabilities
%         pout(k) = P(sir < gamma(k)),
%         where "sir" can represent SIR, SNR or SINR, as described above.
%
% ------------------------ Reference -----------------------------------
% [1] N. C. Beaulieu and C. Cheng, "Efficient Nakagami-m Fading Channel Simulation,"
% IEEE Trans. Veh. Techn., vol. 54, no. 2, pp. 413–424, Mar. 2005.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = SimOutage_BlocksFAMA(Nsamples,gamma, U, mu, L, m, famatype, flagmf, gamma_avg, d0, d, alphaCF, Nant)

if strcmp(famatype, 'Fast') || strcmp(famatype, 'Slow') || strcmp(famatype, 'CFfast') || strcmp(famatype, 'CFslow')

    if strcmp(famatype, 'Fast') || strcmp(famatype, 'Slow')
        Nant = 1;
        d0 = 1;
        d = 1;
        alphaCF = 0;
    end
    mK = m * Nant;
               

    % Pre-allocating
    pout = zeros(size(gamma));

    % Single user flag
    singleUser = (U == 1);

    % Batch size to avoid RAM issues. Nsamples are generated in batches of
    % size batchsize
    batchsize = 1e3;

    % Generate variables by batches
    for batch = 1:batchsize:Nsamples

        % Single user
        if singleUser
            X = GenVariables_SingleUser(batchsize, mu, L, mK);
            sir = max( X ./ ((2*mK)/((1-mu^2)*gamma_avg)) ,[],1);
            
        % Multiuser case
        else
            % Generate variables X and Y
            if  strcmp(famatype, 'Slow')  || strcmp(famatype, 'CFslow')
                [Y, X] = GenVariables(batchsize, m*(U-1), mu, L, mK);
                if gamma_avg == Inf  % sir-based op
                    Y = (d0/d)^alphaCF .* Y;
                else % sinr-based op
                    Y = (d0/d)^alphaCF .* Y + (2*mK)/((1-mu^2)*gamma_avg);
                end

            elseif  strcmp(famatype, 'Fast')  ||  strcmp(famatype, 'CFfast')
                if flagmf % mf = 1, approximated method
                    [Y, X] = GenVariables(batchsize, 1, mu, L, mK);
                    if gamma_avg == Inf  % sir-based op
                        Y = ( (m*(U-1)- m + 1) * (d/d0)^alphaCF  ) .*  Y ;
                    else % sinr-based op
                        beta = (2 * mK  / gamma_avg) * ( ( d / d0 )^alphaCF );
                        msf = ( m * (U-1) + beta + beta^2/(4*m*(U-1)) ) /...
                            ( m * (U-1) - m  + 1 + beta + beta^2/(4*m*(U-1)) );
                        Y = ( (U-1)*(d0/d)^alphaCF + Nant/gamma_avg ) * (m/msf) .* Y;
                    end

                else % mf =~ 1       
                    % for sigma^2 = sigma_s^2 = 1, sqrt of noise power is
                    sigma_eta = sqrt( (2 * mK) / (d0^alphaCF * gamma_avg) );
                    %
                    % Generate Nakagami-m samples, acording with method described in [1]
                    [Y, X] = GenVariablesFast(batchsize, U, m, d, alphaCF, sigma_eta, mu, L, mK);
                    Y = (d0^alphaCF / (1-mu^2)) .* Y;
                end
            end

            % Computes the SIR (or SINR)
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
%   is Nakagami-m distributed.
%   SLOW FAMA and FAST FAMA (for mf = 1)
%---------------------------------------------------------------------
    function [Y, X] = GenVariables(Nsamples, m_slow, mu, L, m)

        % Number of ports
        N = sum(L);

        % Pre-allocating
        X = zeros(N, Nsamples);
        Y = zeros(N, Nsamples);

        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Loop over independent blocks
            for b = 1:length(L)
                b_index = sum(L(1:(b-1)))+1:sum(L(1:b));
                % Correlated samples of channel from/to desired user
                X(b_index,junk_index) = sum((randn(length(b_index),junk_size,m)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,m)).^2 + ...
                    (randn(length(b_index),junk_size,m)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,m)).^2, 3);

                % Correlated samples of channel from/to interferers
                % For slow-FAMA, sum of interferes and fadings (m)
                Y(b_index,junk_index) = sum((randn(length(b_index),junk_size, m_slow)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size, m_slow)).^2 + ...
                    (randn(length(b_index),junk_size, m_slow)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size, m_slow)).^2, 3);
            end

        end

    end % function GenVariables

%---------------------------------------------------------------------
%   Auxiliary function that generates correlated samples for the
%   channel from/to the desired user and the interference. The channel
%   is Nakagami-m.
%   For FAST FAMA
%---------------------------------------------------------------------
    function [Y, X] = GenVariablesFast(Nsamples, U, m, d, alphaCF,sigma_eta, mu, L, mK)

        % Number of ports
        N = sum(L);

        % Pre-allocating
        X = zeros(N, Nsamples);
        Y = zeros(N, Nsamples);

        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Loop over independent blocks
            for b = 1:length(L)
                b_index = sum(L(1:(b-1)))+1:sum(L(1:b));

                % Correlated samples of channel from/to desired user
                X(b_index,junk_index) = sum((randn(length(b_index),junk_size,mK)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,mK)).^2 + ...
                    (randn(length(b_index),junk_size,mK)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,mK)).^2, 3);

                % Correlated samples of channel from interference users
                % Generate Nakagami-m samples, acording with method described in [1]
                for u = 1:((U-1))
                    % One RV to link all ports in a block: randn(1, junk_size)
                    gaux =  sqrt(1-mu^2) * randn(length(b_index), junk_size) + ...
                         mu * randn(1, junk_size) +...
                        1i*(sqrt(1-mu^2)*randn(length(b_index), junk_size) +...
                          mu * randn(1, junk_size));
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
                        % second moment of the Rayleigh RV = 2, then the
                        % CDF of Rayleigh is uniform
                        uniform = ( 1 - exp (- (abs(gaux).^2)./2) );
                        % ancillary variable eta   
                         auxVar = sqrt( log(1./(1-uniform)) ).^(1/m);
                        nkgNormalized = auxVar + ((a1.*auxVar + a2.*(auxVar.^2) + a3.*(auxVar.^3))./(1 + b1.*auxVar + b2.*(auxVar.^2)));
                        % second moment of the Rayleigh RV = 2, then the
                        % denormalized Nakagami-m RV is
                        nkg = sqrt(2) * nkgNormalized ;
                        nkgComplex = nkg .* exp(1i.*angle(gaux));
                    end
                    if m == 1
                        nkgComplex = gaux;
                    end
                    Y(b_index,junk_index) = nkgComplex + Y(b_index,junk_index);
                end

                 % multiply distance 
                if alphaCF ~= 0
                    Y(b_index,junk_index) = d^(-alphaCF/2) * Y(b_index,junk_index);
                end
                % add noise 
                if sigma_eta ~= 0
                   Y(b_index,junk_index) = Y(b_index,junk_index) + (sigma_eta/sqrt(2))...
                      *(randn(length(b_index), junk_size) + 1i*randn(length(b_index), junk_size));
                end
                % absolute value squared
                Y(b_index,junk_index) = abs(Y(b_index,junk_index)).^2;

            end

        end

    end % function GenVariablesFast


%---------------------------------------------------------------------
%   Auxiliary function that generates correlated samples for the
%   channel from/to the desired. 
%   The channel is Nakagami-m distributed.
%   SINGLE USER.
%---------------------------------------------------------------------

    function [X] = GenVariables_SingleUser(Nsamples, mu, L, m)

        % Number of ports
        N = sum(L);

        % Pre-allocating
        X = zeros(N, Nsamples);

        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Loop over independent blocks
            for b = 1:length(L)
                b_index = sum(L(1:(b-1)))+1:sum(L(1:b));
                % Correlated samples of channel from/to desired user
                X(b_index,junk_index) = sum((randn(length(b_index),junk_size,m)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,m)).^2 + ...
                    (randn(length(b_index),junk_size,m)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,m)).^2, 3);

            end

        end

    end % function GenVariables_SingleUser

end

