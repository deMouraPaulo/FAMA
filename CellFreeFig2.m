%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs SIR threshold for a 1D fluid
% antenna in a cell-free FAMA system, with MRT precoding. 
%
% Similar to Fig. 2 of [1], but with d1=d2=d3=500, instead of 
% d1=400, d2=600, d3=800, as in [1].
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


N = 10;                     % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = 3+1;                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation
          
gam = [12,13,14,15,16,17,18];       % SIR threshold (linear scale)
gamdB = pow2db(gam);                % SIR threshold (dB)

m = [1];                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR


cellfreeAlpha = 3;
cellfreeD0 = 200;
cellfreeD = 500;
nAnt = 2;

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_simCFslow = zeros(1,length(gam));
pout_simCFfast = zeros(1,length(gam));
pout_calcCFfast = zeros(1,length(gam));
pout_calcCFslow = zeros(1,length(gam));


% Average correlation coefficient
dm = DeltaMed (N, W);


%-----------------------------------------------------------------
% Outage Probabilities calculation/simulation
%-----------------------------------------------------------------

% User feedback
disp(['Iter ' num2str(1) ' out of ' num2str(2)]);
famatype = "CFslow";
% Calculation
pout_simCFslow(1,:) = SimOutage_BlocksFAMA(Nsamples,gam, U, sqrt(dm), N, m,famatype, Inf, cellfreeD0, cellfreeD, cellfreeAlpha,nAnt);
pout_calcCFslow(1,:) = CalcOutageFAMA(gam, N, dm, U, 'Integral',order, m,famatype, 'SIR', Inf , orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt);

disp(['Iter ' num2str(2) ' out of ' num2str(2)]);
famatype = "CFfast";
pout_simCFfast(1,:) = SimOutage_BlocksFAMA(Nsamples,gam, U, sqrt(dm), N, m,famatype, Inf, cellfreeD0, cellfreeD, cellfreeAlpha,nAnt);
pout_calcCFfast(1,:) = CalcOutageFAMA(gam, N, dm, U, 'Integral',order, m,famatype, 'SIR', Inf , orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt);

 

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

semilogy(gam, pout_simCFslow(1,:), 'b-o', 'linewidth', 2);
hold on; grid on;
semilogy(gam, pout_calcCFslow(1,:), '-v', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2);
semilogy(gam, pout_simCFfast(1,:), '-*', 'linewidth', 2)
semilogy(gam, pout_calcCFfast(1,:), 'r.-','Marker','square', 'linewidth', 2);


leg = legend( "simCFslow"," calcCFslow", "simCFfast", "calcCFfast");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(leg, 'FontSize', 18, 'Interpreter','latex');
xlabel('SIR threshold $\gamma$ ', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-4, 1e-1])
% Titles
str = ['FAMA parameters: $N= $ ',num2str(N),', $W =$ ',...
    num2str(W),', $U = $ ',num2str(U-1), ', $N_{ant} = $ ',num2str(nAnt)];
title(str,'FontSize', 20, 'Interpreter','latex')
str2 = ['$d_0 =$ ',num2str(cellfreeD0),', $d =$ ',num2str(cellfreeD),...
    '$,\alpha$ = ' num2str(cellfreeAlpha), '.  [Fig. 2, ``Cell-free FAMA"]'];
subtitle (str2,'FontSize', 18, 'Interpreter','latex')



