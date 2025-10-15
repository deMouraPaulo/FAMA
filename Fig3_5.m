%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP (SINR-based) vs antenna size (normalized)
% slow FAMA networks under Nakagami-m fading channels.
% 
% Computes: 
% - OP based-SINR: integral expression; 
% - OP based-SINR: Gauss-Laguerre quadrature approximation;
% - OP based-SINR: Monte-Carlo simulation.
%
% Makes Figure 3.5.
% Parameters: m = 2; U = [3, 4]; N = 100; gamma = [2, 3]dB; SNR avg = 20 dB
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

N = 100;                    % Number of ports  
W = 1:10;                    % Antenna size (wavelength normalized)
U = [3, 4];                 % Number of users

min_samples = 5e5;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e3;            % Multiplication factor for simulation 
     
gamdB = [2, 3];              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)

avg_snrdB = 20;                        % SNR, average received (dB)
avg_snr = db2pow( avg_snrdB );        % SNR, average reveiced, linear scale


m = 2;                    % Nakagami-m fading severity
order = 10;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------
pout_calc_int = zeros(length(W), length(gam), length(U));
pout_calc_gl = zeros(length(W), length(gam), length(U));
pout_sim = zeros(length(W), length(gam), length(U));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter = length(W) * length(U) * length(gam);
kiter = 0;

%-----------------------------------------------------------------
% Loop over U
%-----------------------------------------------------------------
for ku = 1:length(U)

    %-----------------------------------------------------------------
    % Loop over W 
    %-----------------------------------------------------------------
    for kw = 1:length(W )

    %-----------------------------------------------------------------
    % Loop over gamma - OP computation
    %-----------------------------------------------------------------
        for kg = 1:length(gam)

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            % Average correlation coefficient
            dm = DeltaMed (N, W(kw));

            % Integral calculation
            pout_calc_int(kw,kg,ku) = CalcOutageFAMA(gam(kg), N, dm, U(ku), 'Integral', order, m, famatype, 'SINR', avg_snr);

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl(kw,kg,ku) = CalcOutageFAMA(gam(kg), N, dm, U(ku), 'Quadrature', order, m,famatype, 'SINR', avg_snr, order);

            % Monte Carlo simulation
            Nsamples = max( round( factor_sim * ( 1 / pout_calc_gl(kw,kg,ku) ) ), min_samples );
            pout_sim(kw,kg,ku) = SimOutage_BlocksFAMA(Nsamples, gam(kg), U(ku), sqrt(dm), N, m,famatype, avg_snr);
        
        end
    end
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

% Colors for each U: "light blue" and "light orange";
colors = [ [71,147,175]; [ 255, 196, 110] ]/255;

%  pout_calc_int(kw,kg,ku)
for kg = length(gam):-1:1
    for ku = length(U):-1:1        
        curveINT(ku) = semilogy(W, squeeze(pout_calc_int(:,kg,ku)), 'Color', colors(ku,:), 'linewidth', 4);
        hold on;
        curveGL = semilogy(W,squeeze(pout_calc_gl(:,kg,ku)), 's', 'Color', 'r', 'linewidth', 4,'MarkerSize', 10);
        curveSIM = semilogy(W, squeeze(pout_sim(:,kg,ku)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    end
end

%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Antenna Size $- \, W$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim([1,10]);
xticks([2,4,6,8,10]);
ylim([0.99e-4, 1]);
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT(1), curveINT(2), curveSIM, curveGL],...
    { ['$U = $ ', num2str(U(1))], ['$U= $ ',num2str(U(2))],'Simulation','Gauss-Laguerre'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='east');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.629,0.594,0.266,0.121],...   
    'String',['Average SNR $=$ ',num2str(avg_snrdB), ' dB', newline ...
    '$N = $ ', num2str(N), newline,...
    '$m =$ ',num2str(m)],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','off');


%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% gamma = 3 dB

% Create arrow
annotation(figure(1),'arrow',[0.3166, 0.36187], [0.735, 0.855],'LineWidth',1, 'HeadLength',15, 'HeadWidth',15)
annotation(figure(1),'arrow',[0.32068,0.34895], [0.6678,0.5548],'LineWidth',1, 'HeadLength',15, 'HeadWidth',15)
% Create textbox
annotation(figure(1),'textbox',...
    [0.2305,0.679,0.117,0.0449],...
    'String',{'$\gamma = 3$ dB'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% gamma = 2 dB

% Create textbox
annotation(figure(1),'textbox',...
    [0.265,0.401,0.133,0.0491],...
    'String',{'$\gamma = 2$ dB'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure(1),'arrow',[0.387,0.471], [0.462,0.7678],'LineWidth',1, 'HeadLength',15, 'HeadWidth',15)
annotation(figure(1),'arrow',[0.366,0.3845],[0.396,0.292],'LineWidth',1, 'HeadLength',15, 'HeadWidth',15)
%------------------------------------------------------------------------



