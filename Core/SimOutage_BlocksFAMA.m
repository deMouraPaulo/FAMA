%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = SimOutage_BlocksNakagami(Nsamples,gamma, U, mu, L, m)
%
% Simulates the outage probability achieved by FAMA under Nakagami-m fading
% with block-diagonal correlation matrix as in Eq. (18). 
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
% - Nsamples: number of Monte-Carlo simulations
% - gamma: SIR threshold (can be a vector)
% - gamma_avg: scalar, average received SNR (used for OP based on SNR and SINR)
% - U: number of users (scalar)
% - mu: correlation factor (scalar) within each block.
% - L: vector containing the size of each block.
% - m: Nakagami-m fading severity
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
%         pout(k) = P(SIR < gamma(k)) according to Eq. (26)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = SimOutage_BlocksFAMA(Nsamples,gamma, U, mu, L, m, famatype, gamma_avg, d0, d, alphaCF, Nant)

if strcmp(famatype, 'Fast') || strcmp(famatype, 'Slow') || strcmp(famatype, 'CFfast') || strcmp(famatype, 'CFslow')


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
            X = GenVariables_SingleUser(batchsize, mu, L, m);
            if strcmp(famatype, 'Slow') || strcmp(famatype, 'Fast')
                sir = max( X ./ ((2*m)/((1-mu^2)*gamma_avg)),[],1);
            elseif strcmp(famatype, 'CFslow') || strcmp(famatype, 'CFfast')
                sir = max( X ./ ((2*m*Nant)/((1-mu^2)*gamma_avg)) ,[],1);
            end
        % Multiuser case
        else
            % Generate variables X and Y according to Eq. (26)
            if strcmp(famatype, 'Slow')
                [Y, X] = GenVariables(batchsize, m*(U-1), mu, L, m);
                Y = Y + (2*m)/((1-mu^2)*gamma_avg);
            elseif strcmp(famatype, 'Fast')
                [Y, X] = GenVariables(batchsize, 1, mu, L, m);
                Y = ( m*(U-2)+1) .*  Y + (2*m)/((1-mu^2)*gamma_avg);
            elseif strcmp(famatype, 'CFslow')
                [Y, X] = GenVariables(batchsize, m*(U-1), mu, L, m*Nant);
                Y = (d0/d)^alphaCF .* Y + (2*m*Nant)/((1-mu^2)*gamma_avg);
            elseif strcmp(famatype, 'CFfast')
                [Y, X] = GenVariables(batchsize, 1, mu, L, m*Nant);
                Y =  (m*(U-2)+1)*(d0/d)^alphaCF .* Y + (2*m*Nant)/((1-mu^2)*gamma_avg);
            end

            % Computes the SIR according to Eq. (26)
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
%   is Nakagami-m distributed with correlation matrix given by
%   Sigma = V'*diag(lambda)*V
%   Based on s-FAMA, but can be used for other types of FAMA
%---------------------------------------------------------------------
    function [Y, X] = GenVariables(Nsamples, Util, mu, L, m)

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
                Y(b_index,junk_index) = sum((randn(length(b_index),junk_size, Util)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size, Util)).^2 + ...
                    (randn(length(b_index),junk_size, Util)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size, Util)).^2, 3);
            end

        end

    end % function GenVariables


%---------------------------------------------------------------------
%   Auxiliary function that generates correlated samples for the
%   channel from/to the desired. 
%   The channel is Nakagami-m distributed with correlation matrix given by
%   Sigma = V'*diag(lambda)*V
%   Based on s-FAMA, but can be used for other types of FAMA
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

