%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP (SINR-based) vs Number of ports for 
% slow FAMA networks under Nakagami-m fading channels.
% 
% Computes: 
% - OP based-SINR: integral expression; 
% - OP based-SINR: Gauss-Laguerre quadrature approximation;
% - OP based-SINR: Monte-Carlo simulation.
%
% Makes Figure 3.4.
% Parameters: m = [1,2,3]; U = [3, 5]; W = 1; gamma = 1 dB; SNR avg = 20 dB
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

famatype = "Slow";     %  FAMA type

N= round(logspace(1, log10(150), 10 )); % Number of ports  
Nsimul = N;
W = 1;                      % Antenna size (wavelength normalized)
U = [3, 5];                 % Number of users

min_samples = 5e5;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e2;            % Multiplication factor for simulation 
     
gamdB = 1;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)

avg_snrdB = 20;                        % SNR, average received (dB)
avg_snr = db2pow( avg_snrdB );        % SNR, average reveiced, linear scale


m = [1,2,3];                    % Nakagami-m fading severity
order = 10;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------
pout_calc_int = zeros(length(N), length(m), length(U));
pout_calc_gl = zeros(length(Nsimul), length(m), length(U));
pout_sim = zeros(length(Nsimul), length(m), length(U));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter = length(N) *  length(m) * length(U);
kiter = 0;

%-----------------------------------------------------------------
% Loop over  m
%-----------------------------------------------------------------
for km = 1:length(m)

    %-----------------------------------------------------------------
    % Loop over U
    %-----------------------------------------------------------------
    for ku = 1:length(U)

        %-----------------------------------------------------------------
        % Loop over N - OP computation
        %-----------------------------------------------------------------
        for kn = 1:length(N )

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);


            % Average correlation coefficient
            dm = DeltaMed (N(kn), W);

            % Integral calculation
            pout_calc_int(kn,km,ku) = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Integral', order, m(km), famatype, 'SINR', avg_snr);        
            
            
        end

        %-----------------------------------------------------------------
        % Loop over Nsimul - Simulation and Quadrature Calculation
        %-----------------------------------------------------------------
        for kn = 1:length( Nsimul )


            % Average correlation coefficient
            dm = DeltaMed (Nsimul(kn), W);

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl(kn,km,ku) = CalcOutageFAMA(gam, Nsimul(kn), dm, U(ku), 'Quadrature', order, m(km),famatype, 'SINR', avg_snr, order);

            % Monte Carlo simulation
            Nsamples = max( round( factor_sim * ( 1 / pout_calc_gl(kn,km,ku) ) ), min_samples );
            pout_sim(kn,km,ku) = SimOutage_BlocksFAMA(Nsamples, gam, U(ku), sqrt(dm), Nsimul(kn), m(km),famatype, avg_snr);

            
        end
    end 
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

% Colors for each "m" curve: "light blue", "light organge and "red"";
colors = [ [71,147,175]; [ 255, 196, 110];  [255, 0, 0] ]/255;

for km = length(m):-1:1
    for ku = length(U):-1:1
        curveINT(km) = loglog(N, squeeze(pout_calc_int(:,km,ku)), 'Color', colors(km,:), 'linewidth', 4);
        hold on;
        curveGL = loglog(Nsimul,squeeze(pout_calc_gl(:,km,ku)), 's', 'Color', 'r', 'linewidth', 4,'MarkerSize', 10);
        curveSIM = loglog(Nsimul, squeeze(pout_sim(:,km,ku)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    end
end

%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Ports $- \, N$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [10, 150] );
ylim( [0.99e-5, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT(1), curveINT(2), curveINT(3), curveGL, curveSIM],...
    { ['$m = $ ', num2str(m(1))], ['$m= $ ',num2str(m(2))], ['$m= $ ',num2str(m(3))],'Gauss-Laguerre','Simulation'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='west');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.464,0.133,0.3161,0.1471],...   
    'String',['Average SNR $=$ ',num2str(avg_snrdB), ' dB', newline ...
    '$\gamma = $ ', num2str(gamdB), ' dB' newline,...
    '$W =$ ',num2str(W)],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'BackgroundColor','w',...
    'FontSize',24,...
    'FitBoxToText','off');


%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% U = 5
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.7983,0.720,0.0394,0.205]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.783,0.667,0.0826,0.0447],...
    'String',{'$U = 5$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% U = 3
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.540,0.460,0.0394,0.316]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.511,0.405,0.0826,0.0447],...
    'String',{'$U = 3$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------


