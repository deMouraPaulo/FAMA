%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots Outage probability vs Number of ports for a network with:
% - cell-free with maximum ratio transmission (MRT) precoding 
% at base station and users with fluid antennas system (FAS);
% - fast and slow fluid antenna multiple access (FAMA)
% under Nakagami-m fading channels;
% - single fixed position antenna, for comparative analysis;
% - constant correlation model.
% 
% Computes: 
% - OP based on SIR using constant correlation model. 
% 
% Makes Figure 5.2.
%
% Parameters:
% U = 4; nAnt = [2,4]; nAnt_fix(compare with FPA) = [15,20]; 
% m = 2; gamdB = 0; W = 5; d0 = 100; d = 100; alphaCF = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------
clc
close all
clear

%addpath('Core/')
addpath(fullfile('..','Core'))


%-------------------------------------------------------------------------
% Flags
%-------------------------------------------------------------------------

% set flags for f-FAMA calculus (for s-Fama don't care)

% flag, approx. mf = 1
flagmf = false;

% flag, 2nd integral
flag2int = false;

%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------


N = [1,2,3,4,5,7,8,9,10,12,15,20];      % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = 3+1;                    % Number of interference BS + 1

gam = 1;                    % SIR threshold (linear scale)
gamdB = pow2db(gam);        % SIR threshold (dB)



m = 2;                    % Nakagami-m fading serverity
order = 30;                 % Orderm of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR


cellfreeAlpha = 3;
cellfreeD0 = 100;
cellfreeD = 100;

% Number of BS antennas, when FAS is used in the reception.
nAnt = [2,4];

% Number of BS antennas, when one fixed antenna is used in the reception.
nAnt_fix = [15,20];

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_CFslow = zeros(length(N),length(nAnt));
pout_CFfast = zeros(length(N),length(nAnt));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%-----------------------------------------------------------------

   
% Loop number of Antenas (nAnt)
for ka = 1:length(nAnt)

    % User feedback
    disp(['Iter ' num2str(ka) ' out of ' num2str(length(nAnt))]);

    % Loop number of ports (N)
    for kn = 1:length(N )                
        if N(kn) == 1
            pout_CFfast(kn,ka) = pout1ant_fast(m,nAnt(ka),U,gam,cellfreeD0,cellfreeD,cellfreeAlpha);
            pout_CFslow(kn,ka) = pout1ant_slow(m,nAnt(ka),U,gam,cellfreeD0,cellfreeD,cellfreeAlpha);
        else
            % Constant correlation coefficient
            dm = DeltaMed (N(kn), W);
            famatype = "CFfast";
            pout_CFfast(kn,ka) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral',order, m,famatype,flag2int,flagmf,Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt(ka));
            famatype = "CFslow";
            pout_CFslow(kn,ka) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral',order, m,famatype, flag2int,flagmf, Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt(ka));
        end
    end
end

% One fixed antenna in the reception
pout_fix = zeros(length(N),length(nAnt_fix));
for ka = 1:length(nAnt_fix)  
    famatype = "CFslow";
    pout_fix(:,ka) = pout1ant_slow(m,nAnt_fix(ka),U,gam,cellfreeD0,cellfreeD,cellfreeAlpha);
end



%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

Fig5_2curves


% -------------------- Single antenna SLOW -----------------
function [pout] = pout1ant_slow (m,K,U,gam,d0,d,alf)
    U = U-1;
    fun = @(r) r.^(m*U-1) .* exp(-r) .* gammainc( (d0./d).^alf .* gam.* r, m.*K);
    pout = integral(fun,0,inf)./(gamma(m.*U));
end


% -------------------- Single antenna FAST -----------------
function [pout] = pout1ant_fast (m,K,U,gam,d0,d,alf)
    U = U-1;
    fun = @(r)   exp(-r) .* gammainc((d0./d).^alf .* (m .* U - m + 1).* gam .* r, m.*K);
    pout = integral(fun,0,inf);
end