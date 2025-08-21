%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = CalcOutageNakagami(gamma_v, L, rho, U, method, order, m)
%
% Calculates upper bound OP of fast and slow FAMA for 
% the block-diagonal correlation approximation 
% 
% Parameters:
%
% - gamma_v: vector containing the SIR thresholds in linear scale. If U ==
%            1, then gamma_v represents the SNR threshold.
% - B: number of blocks, or eigenvalues (scalar)
% - U: number of users (scalar). If U ==1, then Eq. (43) is computed
% - m: fading
% - pout: vector containing the OP (same size as gamma_v), i.e., 
%         pout(k) = P(SIR < gamma_v(k))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [pout] = OutageUBblocks(gamma_v, B, U, mnkg, famatype)

    if strcmp(famatype, 'Fast')
        % for fast-FAMA = mtil -> Util 
        mtil =  1 / ( 1 - (mnkg-1) / (mnkg * (U-1)) );
        Util = mtil;
        % for slow-FAMA, gamma_th -> gamma_th * Uhat
        Uhat = mnkg * (U - 2) + 1;
        gamma_v =  Uhat .* gamma_v;
    
    elseif strcmp(famatype, 'Slow')
        Util = mnkg * (U-1);
        
    else
         error('Type of FAMA must be indicaded: Fast or Slow')
    end

    sumatory = zeros(1,length(gamma_v));
    for i = 0:(mnkg-1)
        aux = (gamma_v.^i .* pochhammer(Util, i) )./...
            (factorial(i).*(gamma_v + 1).^(Util + i)  );
        sumatory = aux + sumatory(1,:);
    end
    pout = (1 - sumatory).^B;

end
