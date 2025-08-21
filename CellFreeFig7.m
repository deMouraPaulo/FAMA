%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs num transm. antennas for a 1D fluid
% antenna in a cell-free FAMA system, with MRT precoding. 
%
% Similar to Fig. 7 of [1].
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


N = [25,50];                     % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = [2,3,4,5];                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

     
gam = 3;                % SIR threshold (linear scale)
gamdB = db2pow(gam);               % SIR threshold (dB)


m = [1];                      % Nakagami-m fading serverity
order = 30;                    % Order of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR

cellfreeAlpha = 3;
cellfreeD0 = 100;
cellfreeD = 100;
nAnt = 2;

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_CFfast = zeros(length(U),length(N));
pout_CFslow = zeros(length(U),length(N));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

% Loop over  N
for kn = 1:length(N)

    % User feedback
    disp(['Iter ' num2str(kn) ' out of ' num2str(length(N))]);

    % Average correlation coefficient
    dm = DeltaMed (N(kn), W);

    % Loop over interf U
    for ku = 1:length(U )

        famatype = "CFfast";
        % CalcOutageFAMA(gamma_v, L, rho, U, method, order, m, famatype,opbasedSINR, gamma_avg, orderSINR, d0, d, alphaCF, Nant)
        pout_CFfast(ku,kn) = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt);
        famatype = "CFslow";
        pout_CFslow(ku,kn) = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Integral',order, m,famatype, 'SIR', Inf, orderSINR, cellfreeD0, cellfreeD, cellfreeAlpha, nAnt);

    end
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for kn = length(N):-1:1

semilogy(U-1, pout_CFfast(:,kn), '-vb',  'linewidth', 2);
hold on; grid on;
semilogy(U-1, pout_CFslow(:,kn), '-*r', 'linewidth', 2);

end

l = legend( "CF-fast", "CF-slow");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('Number of interfering BS ', 'FontSize', 20, 'Interpreter','latex');
xticks(U-1);
ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-20, 1e0])
str = [' FAMA parameters: $W = $ ',num2str(W), ', $\gamma = $ ',num2str(gam) ];
title(str,'FontSize', 20, 'Interpreter','latex')
str2 = ['$d_0 =$ ',num2str(cellfreeD0),', $d =$ ',num2str(cellfreeD),...
    '$,\alpha$ = ' num2str(cellfreeAlpha),'.  [Fig. 7, ``Cell-free FAMA"]'];
subtitle (str2,'FontSize', 18, 'Interpreter','latex')



