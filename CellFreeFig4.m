%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs num transm. antennas for a 1D fluid
% antenna in a cell-free FAMA system, with MRT precoding. 
%
% Similar to Fig. 4 of [1]
% [1] T. Han, Y. Zhu, K.-K. Wong, G. Zheng, and H. Shin, 
% "Cell-free Fluid Antenna Multiple Access Networks," 
% IEEE Trans. Wireless Commun., pp. 1–1, Apr. 2025,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------
clc
close all
clear

addpath('Core/')

%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------


N = [1,2,5,7,10,15,20];                     % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = 3+1;                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

%gam = linspace(12,18, 7);     % SIR threshold (dB)
% gam = 10.^(gamdB/10);            % SIR threshold (linear scale)
gam = [1];
gamdB = pow2db(gam);

m = [1];                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR


cellfreeAlpha = 3;
cellfreeD0 = 100;
cellfreeD = 100;
nAnt = [2,4];

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
        
        % Average correlation coefficient 
        if N(kn) ~= 1
            dm = DeltaMed (N(kn), W);
        elseif N(kn) == 1
            dm = 0.99999;
        end
        famatype = "CFfast";
        pout_CFfast(kn,ka) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt(ka));
        famatype = "CFslow";
        pout_CFslow(kn,ka) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt(ka));

    end
end


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for ka = length(nAnt):-1:1
    loglog(N, pout_CFfast(:,ka), '-vb',  'linewidth', 2);
    hold on; grid on;
    loglog(N, pout_CFslow(:,ka), '-*r', 'linewidth', 2);
end


l = legend( "CFfast", "CFslow");
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');

xlabel('Number of ports ', 'FontSize', 20, 'Interpreter','latex');
xticks(N);

ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-10, 1e0])

str = [' FAMA parameters:  $W =$ ', num2str(W),...
    ', $U_{\mathrm{interf}} = $ ',num2str(U-1), ', $\gamma = $',num2str(gam)];
title(str,'FontSize', 20, 'Interpreter','latex')
str2 = ['$d_0 =$ ',num2str(cellfreeD0),', $d =$ ',num2str(cellfreeD),...
    '$,\alpha$ = ' num2str(cellfreeAlpha),'.  [Fig. 4, ``Cell-free FAMA"]'];
subtitle (str2,'FontSize', 18, 'Interpreter','latex')




