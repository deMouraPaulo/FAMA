%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function mg = MultiplexGain(num_pool, num_users, out_prob ) 
%
% DESCRIPTION:
%   Calculates the multiplexing gain based on the exact expression
%   expressed in eq. (21a) (paper opportunistic FAMA).
%
% INPUT:
%   num_pool - (int) Number of users in the pool (M)
%   num_users - (int) Number of users selected (U)
%   out_prob - (double) Outage Probability (p)
%
% OUTPUT:
%   mg - (double)  multiplexing gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mg = MultiplexGain(num_pool, num_users, out_prob )

mg = 0;
for u = (num_pool - num_users + 1):num_pool
  mg = mg + betainc(1-out_prob, num_pool - u + 1, u);
end