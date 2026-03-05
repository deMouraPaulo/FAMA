%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs Average SNR for 
% slow FAMA networks under Nakagami-m fading channels.
%
% Computes: 
% - OP based-SIR: integral expression; 
% - OP based-SINR: integral expression; 
% - OP based-SINR: Gauss-Laguerre quadrature approximation;
% - OP based-SINR: Monte-Carlo simulation.
% 
% Makes Figure 3.2.
% Parameters: m = [1,2,3]; N = 150; W = 1; U = 4; gamma = 1 dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

clc
clear
close all
addpath('Core/')


%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------

famatype = 'Slow';     % Set FAMA type

N = 150;               % Number of ports
W = 1;                 % Antenna size (wavelength normalized)
U = 4;                 % Number of users

Nsamples = 5e5;        % Number of samples for Monte-Carlo simulation
     
gamdB = 1;              % SINR threshold (dB)
gam = db2pow(gamdB);    % SINR threshold (linear scale)

avg_snrdB = linspace( -10, 25, 20 );  % SNR, average received (dB)
avg_snr = db2pow( avg_snrdB );        % SNR, average reveiced, linear scale


m = [1,2,3];              % Nakagami-m fading severity
order = 30;               % Order of GL quadrature
orderSINR = 30;           % Order of 1st GL quadrature for OP-SINR


%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_sir = zeros(length(m),length(avg_snr));
pout_calc_int = zeros(length(m),length(avg_snr));
pout_calc_gl = zeros(length(m),length(avg_snr));
pout_sim = zeros(length(m),length(avg_snr));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

% Average correlation coefficient
dm = DeltaMed (N, W);

% Loop over  m
for km = 1:length(m)


    % User feedback
    disp(['Iter "m": ' num2str(km) ' out of ' num2str(length(m))]);

    % Loop over avg_snr
    for kg = 1:length(avg_snr )

        % User feedback
        disp(['Iter "SNR": ' num2str(kg) ' out of ' num2str(length(avg_snr))]);

        pout_calc_gl(km,kg) = CalcOutageFAMA(gam, N, dm, U, 'Quadrature', order, m(km),famatype, 'SINR', avg_snr(kg), orderSINR);
        pout_calc_int(km,kg) = CalcOutageFAMA(gam, N, dm, U, 'Integral', order, m(km),famatype, 'SINR', avg_snr(kg), orderSINR);

        % Monte Carlo simulation
        NsamplesAdjusted = max( round( 1000 * ( 1 / pout_calc_int(km,kg) ) ), Nsamples );
        pout_sim(km,kg) = SimOutage_BlocksFAMA(NsamplesAdjusted,gam, U, sqrt(dm), N, m(km),famatype, avg_snr(kg));

    end
    pout_calc_sir(km,:) = CalcOutageFAMA(gam, N, dm, U, 'Quadrature', order, m(km),famatype, 'SIR', Inf, orderSINR);

end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for km = length(m):-1:1

    curveSIR = semilogy(avg_snrdB, pout_calc_sir(km,:), '--', 'Color','r', 'Linewidth', 3);
    hold on; 
    curveINT = semilogy(avg_snrdB, pout_calc_int(km,:), 'Color', [71,147,175]/255, 'linewidth', 4);
    curveGL = semilogy(avg_snrdB, pout_calc_gl(km,:), 's', 'Color', 'r', 'linewidth', 2, 'MarkerSize', 10);
    curveSIM = semilogy(avg_snrdB, pout_sim(km,:), '+', 'Color', 'k', 'linewidth', 2,'MarkerSize',8);

end

%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Average SNR (dB)', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [-10, 25] );
ylim( [1e-3, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

% Create Legend
legend( [curveSIR, curveINT, curveSIM, curveGL],...
    {'SIR', 'SINR', 'Simulation', 'Gauss-Laguerre'}, ...
    'Interpreter', 'Latex', 'FontSize', 24, Location='east' );

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.144,0.396,0.129,0.1778],...
    'String',['$N = $ ', num2str(N), newline,...
    '$W=$ ',num2str(W), newline,...
    '$U=$ ',num2str(U),newline,...
    '$\gamma= $ ',num2str(gamdB),' dB'],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','on');

% Create ellipse m = 1
annotation(figure(1),'ellipse',...
    [0.5236, 0.216, 0.0264, 0.135]);
% Create ellipse m = 2
annotation(figure(1),'ellipse',...
    [0.566, 0.579, 0.0240, 0.0937]);
% Create ellipse m = 3
annotation(figure(1),'ellipse',...
    [0.650, 0.745, 0.0200, 0.0724]);

% Create textbox m = 1
annotation(figure(1),'textbox',...
    [0.438, 0.240, 0.1, 0.0515],...
    'String',{'$m=1$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'EdgeColor',[1 1 1]);

% Create textbox m = 2
annotation(figure(1),'textbox',...
    [0.477, 0.5469, 0.1 0.0515],...
    'String',{'$m=2$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'EdgeColor',[1 1 1]);

% Create textbox m = 3
annotation(figure(1),'textbox',...,
    [0.621, 0.810 0.1 0.0515],...
    'String',{'$m=3$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'EdgeColor',[1 1 1]);
