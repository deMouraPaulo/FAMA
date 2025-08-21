%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs num transm. antennas for a 1D fluid
% antenna in a cell-free FAMA system, with MRT precoding. 
%
% Similar to Fig. 3 of [1], but with d1=d2=d3=500, instead of 
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

N = 2;                     % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = 3+1;                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

gam = [14,18];              % SIR threshold (linear scale)
gamdB = pow2db(gam);         % SIR threshold (dB)

m = [1];                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR


cellfreeAlpha = 3;
cellfreeD0 = 200;
cellfreeD = 500;
nAnt =  1:1:10;

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_CFfast = zeros(length(nAnt),length(gam));
pout_CFslow = zeros(length(nAnt),length(gam));

% Average correlation coefficient
dm = DeltaMed (N, W);


%-----------------------------------------------------------------
% Outage Probabilities calculation
%-----------------------------------------------------------------

% User feedback
disp(['Iter ' num2str(1) ' out of ' num2str(2)]);
famatype = "CFfast";

for kg = 1:length(gam)
    for ka = 1:length(nAnt)
        pout_CFfast(ka,kg) = CalcOutageFAMA(gam(kg), N, dm, U, 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt(ka));
    end
end

disp(['Iter ' num2str(2) ' out of ' num2str(2)]);
famatype = "CFslow";
for kg = 1:length(gam)
    for ka = 1:length(nAnt)
        pout_CFslow(ka,kg) = CalcOutageFAMA(gam(kg), N, dm, U, 'Integral',order, m,famatype,'SIR', Inf, orderSINR,cellfreeD0, cellfreeD, cellfreeAlpha, nAnt(ka));
    end
end


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for kg = length(gam):-1:1
    semilogy(nAnt, pout_CFfast(:,kg), 'b--o', 'linewidth', 2);
    hold on; grid on;
    semilogy(nAnt, pout_CFslow(:,kg), 'r.-','Marker','square', 'linewidth', 2);
end

leg = legend( "CF-fast", "CF-slow");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(leg, 'FontSize', 18, 'Interpreter','latex');
xlabel('Number of antennas in a BS', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-4, 1])
%
% Titles
str = [' FAMA parameters: $N= $ ',num2str(N),', $W =$ ',num2str(W),', $U = $ ',num2str(U)];
title(str,'FontSize', 20, 'Interpreter','latex')
str2 = ['$d_0 =$ ',num2str(cellfreeD0),', $d =$ ',num2str(cellfreeD),...
    '$,\alpha$ = ' num2str(cellfreeAlpha),'.  [Fig. 3, ``Cell-free FAMA"]'];
subtitle (str2,'FontSize', 18, 'Interpreter','latex')



