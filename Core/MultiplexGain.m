%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function mg = MultiplexGain(num_pool, num_users, out_prob ) 
%
% DESCRIPTION:
%   Calculates the multiplexing gain based on the exact expression
%   expressed in Eq. (21a) [*].
%
% INPUT:
%   num_pool - (int) Number of users in the pool (M)
%   num_users - (int) Number of users selected (U)
%   out_prob - (double) Outage Probability (p)
%
% OUTPUT:
%   mg - (double)  multiplexing gain
%
% [*] K. -K. Wong, K. -F. Tong, Y. Chen, Y. Zhang and C. -B. Chae, 
% "Opportunistic Fluid Antenna Multiple Access", 
% IEEE Trans. Wireless Commun., vol. 22, no. 11, pp. 7819-7833, 
% Nov. 2023, doi: 10.1109/TWC.2023.3255940.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mg = MultiplexGain(num_pool, num_users, out_prob )

mg = 0;
for u = (num_pool - num_users + 1):num_pool
  mg = mg + betainc(1-out_prob, num_pool - u + 1, u);
end