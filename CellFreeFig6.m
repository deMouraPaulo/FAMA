%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs num transm. antennas for a 1D fluid
% antenna in a cell-free FAMA system, with MRT precoding. 
%
% Similar to Fig. 6 of [1]
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


N = 60;                     % Number of ports

W = [1,2,3,4,5];                       % Antenna size (wavelength normalized)
U = [2,3];                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

     
gam = 6;                % SIR threshold (linear scale)
gamdB = db2pow(gam);               % SIR threshold (dB)



m = [1];                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR


cellfreeAlpha = 3;
cellfreeD0 = 100;
cellfreeD = 100;
nAnt = 2;

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_CFfast = zeros(length(U),length(W));
pout_CFslow = zeros(length(U),length(W));

%-----------------------------------------------------------------
% Outage Probabilities calculation
%-----------------------------------------------------------------

% Loop over U
for ku = 1:length(U)

% User feedback
disp(['Iter ' num2str(ku) ' out of ' num2str(length(U))]);


% Loop over  W
for kw = 1:length(W )

   
    % Average correlation 
    dm = DeltaMed (N, W(kw));
    famatype = "CFfast";
    pout_CFfast(ku,kw) = CalcOutageFAMA(gam, N, dm, U(ku), 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt);
     famatype = "CFslow";
    pout_CFslow(ku,kw) = CalcOutageFAMA(gam, N, dm, U(ku), 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt);

end
end
%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for ku = length(U):-1:1

    semilogy(W, pout_CFfast(ku,:), '-vb',  'linewidth', 2);
    hold on; grid on;
    semilogy(W, pout_CFslow(ku,:), '-*r', 'linewidth', 2);

end

l = legend( "CF-fast", "CF-slow");
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
%
xlabel('W (in wavelenghs) ', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-8, 1e0])
%
str = [' FAMA parameters:  $N =$ ', num2str(N),...
    ', $\gamma = $ ',num2str(gam) ,', $N_{ant} = $ ',num2str(nAnt)];
title(str,'FontSize', 20, 'Interpreter','latex')
str2 = ['$d_0 =$ ',num2str(cellfreeD0),', $d =$ ',num2str(cellfreeD),...
    '$,\alpha$ = ' num2str(cellfreeAlpha),'.  [Fig. 6, ``Cell-free FAMA"]'];
subtitle (str2,'FontSize', 18, 'Interpreter','latex')


