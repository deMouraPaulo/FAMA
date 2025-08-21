% Absolute average of the correlation coefficients for FAS, as in
% [Ref.1, Eq.4].
% Obs.: dm = mu^2 [1, Eq.4].
%
function [dm] = DeltaMed (n,w)
dm=0; 
for k = 1:(n-1)
    dm = dm + (n-k) * besselj( 0, ( (2*pi*k/(n-1)) * w ) );
end
dm = abs( dm * (2/(n*(n-1)) ) );
end
%
%
% Reference
% [Ref. 1] K. -K. Wong, K. F. Tong, Y. Chen, and Y. Zhang, "Closed-form
% expressions for spatial correlation parameters for performance
% analysis of fluid antenna systems,” Electron. Lett., vol. 58, no.11, pp. 454–457, May 2022.
