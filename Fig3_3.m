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
% Makes Figure 3.3.
% Parameters: m = 2; U = [3,4,5]; N = [100,150]; W = 1; gamma = 0 dB.
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

N = [100,150];               % Number of ports
W = 1;                 % Antenna size (wavelength normalized)
U = [3,4,5];                 % Number of users

Nsamples = 5e5;        % Number of samples for Monte-Carlo simulation
     
gamdB = 0;              % SINR threshold (dB)
gam = db2pow(gamdB);    % SINR threshold (linear scale)

avg_snrdB = linspace( -10, 25, 30 );  % SNR, average received (dB)
avg_snr = db2pow( avg_snrdB );        % SNR, average reveiced, linear scale

avg_snr_simdB =  linspace( -10, 25, 5 );
avg_snr_sim = db2pow( avg_snr_simdB );



m = 2;                    % Nakagami-m fading severity
order = 30;               % Order of GL quadrature
orderSINR = 30;           % Order of 1st GL quadrature for OP-SINR


%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_sir = zeros(length(N), length(U),length(avg_snr));
pout_calc_int = zeros(length(N), length(U),length(avg_snr));
pout_calc_gl = zeros(length(N), length(U),length(avg_snr_sim));
pout_sim = zeros(length(N), length(U),length(avg_snr_sim));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter = length(N) * length(U) * length(avg_snr);
kiter = 0;
% Loop over  N
for kn = 1:length(N)

    % User feedback
    disp(['Iter "N": ' num2str(kn) ' out of ' num2str(length(N))]);

    % Average correlation coefficient
    dm = DeltaMed (N(kn), W);

    % Loop over  U
    for ku = 1:length(U)

        % User feedback
        disp(['Iter "U": ' num2str(ku) ' out of ' num2str(length(U))]);

        % OP-based SIR
        pout_calc_sir(kn,ku,:) = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Quadrature', order, m,famatype, 'SIR', Inf, orderSINR);

        % Loop over avg_snr
        for kg = 1:length(avg_snr )

            % User feedback
            disp(['Iter "SNR": ' num2str(kg) ' out of ' num2str(length(avg_snr))]);
            kiter = kiter + 1;
            disp(['Total Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            % Integral calculation
            pout_calc_int(kn,ku,kg) = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Integral', order, m,famatype, 'SINR', avg_snr(kg), orderSINR);

        end

        % Loop over snr_sim
        for kg = 1:length(avg_snr_sim )
            
            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl(kn,ku,kg) = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Quadrature', order, m,famatype, 'SINR', avg_snr_sim(kg), orderSINR);
            
            % Monte Carlo simulation
            NsamplesAdjusted = max( round( 100 * ( 1 / pout_calc_gl(kn,ku,kg) ) ), Nsamples );
            pout_sim(kn,ku,kg) = SimOutage_BlocksFAMA(NsamplesAdjusted,gam, U(ku), sqrt(dm), N(kn), m,famatype, avg_snr_sim(kg));
        
        end

    end

end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)


% Colors for each "N" curve: "light organge and "light blue",;
colors = [  [ 255, 196, 110];  [71,147,175]; ]/255;

for kn = length(N):-1:1
    for ku = length(U):-1:1
        curveSIR = semilogy(avg_snrdB, squeeze(pout_calc_sir(kn,ku,:)), '--', 'Color','r', 'Linewidth', 3);
        hold on;
        curveINT(kn) = semilogy(avg_snrdB, squeeze(pout_calc_int(kn,ku,:)), 'Color', colors(kn,:), 'linewidth', 4);
        curveGL = semilogy(avg_snr_simdB,squeeze(pout_calc_gl(kn,ku,:)), 's', 'Color', 'r', 'linewidth', 2, 'MarkerSize', 10);
        curveSIM = semilogy(avg_snr_simdB, squeeze(pout_sim(kn,ku,:)), '+', 'Color', 'k', 'linewidth', 2,'MarkerSize', 10);
    end
end

%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Average SNR (dB)', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [-10, 25] );
ylim( [1e-6, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveSIR, curveINT(1), curveINT(2), curveGL, curveSIM],...
    {'SIR', ['SINR, $N = $ ', num2str(N(1))],...
    ['SINR, $N = $ ',num2str(N(2))], 'Gauss-Laguerre','Simulation'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',2, Location='south');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.1576,0.4896,0.1179,0.139],...   
    'String',['$m = $ ', num2str(m), newline,...
    '$W=$ ',num2str(W), newline,...
    '$\gamma= $ ',num2str(gamdB),' dB'],...
    'Interpreter','latex',...
    'FontSize',24,...
    'HorizontalAlignment','center',...
    'FitBoxToText','on');

%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

% Create ellipse U = 3
annotation(figure(1),'ellipse',...
    [0.74, 0.325, 0.028, 0.15,]);
% Create ellipse U = 4
annotation(figure(1),'ellipse',...
    [0.74, 0.642, 0.024, 0.0969]);
% Create ellipse U = 5
annotation(figure(1),'ellipse',...
    [0.74, 0.805, 0.0208, 0.0714]);


% Create textbox U = 3
annotation(figure(1),'textbox',...
    [0.78, 0.455, 0.0846, 0.0373],...
    'String','$U = 3$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox U = 4
annotation(figure(1),'textbox',...
    [0.78 0.735, 0.0846, 0.0373],...
    'String','$U = 4$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox U = 5
annotation(figure(1),'textbox',...
    [0.78, 0.867, 0.0846, 0.0373],...
    'String',{'$U = 5$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

