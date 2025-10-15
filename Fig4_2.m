%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs SIR threshold gamma(dB)
% Fast and slow FAMA networks with block-correlation model.
% 
% Makes Figure 4.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

tic
clc
clear
close all
addpath('Core/')


%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------

famatype = ["Slow";"Fast"];     %  FAMA type
N = 100;                   % Number of ports
W = 1;                      % Antenna size (wavelength normalized)
U = [5, 10];                % Number of users

min_samples = 1e4;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e2;            % Multiplication factor for simulation 
     
gamdB = -15:10;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)

gamdB_sim =  linspace(-15,10,10);   % SIR threshold (dB) for simulation
gam_sim = db2pow(gamdB_sim);        % SIR threshold (linear scale) for simulation


m = [1,2,3];                    % Nakagami-m fading severity
order = 50;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_int_blocks = zeros(length(gam ),length(U), length(m), length(famatype));
pout_calc_gl_blocks = zeros(length(gam_sim ),length(U), length(m), length(famatype));
pout_sim_blocks = zeros(length(gam_sim ),length(U), length(m), length(famatype));


%-----------------------------------------------------------------
%  Block correlation
%----------------------------------------------------------------
% Correlation matriz
Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));
% Eigenvalues
rho = sort(eig(Sigma_jakes),'descend');
% Correlation coefficent per block
deltab = 0.97;
% Number of domminant eigenvalues
Num_eig = sum(rho > N/100);

% Algorithm 1. L: vector with block sizes (Lb)
L = BlockCorrelation(N, rho, Num_eig, deltab);


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter =  length(U)  *  length(m) * length(famatype);
kiter = 0;

%-----------------------------------------------------------------
% Loop over  m
%-----------------------------------------------------------------
for km = 1:length(m)

    %-----------------------------------------------------------------
    % Loop over  U
    %-----------------------------------------------------------------
    for ku = 1:length(U)

        %-----------------------------------------------------------------
        % Loop over famatype 
        %-----------------------------------------------------------------
        for kfama = 1:length(famatype )

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);
     
            % Integral calculation
             pout_calc_int_blocks(:,ku,km,kfama) = CalcOutageFAMA(gam, L, deltab, U(ku), 'Integral', order, m(km), famatype(kfama), 'SIR', Inf);

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl_blocks(:,ku,km,kfama) = CalcOutageFAMA(gam_sim, L, deltab, U(ku), 'Quadrature', order, m(km),famatype(kfama), 'SIR', Inf);

            %-------------------------------------------------------------
            % Monte Carlo simulation
            % Loop over gam_sim
            %-------------------------------------------------------------

            for kg = 1:length (gam_sim)
                % Pout calculated
                pout = pout_calc_gl_blocks(kg,ku,km,kfama);

                % Avoid simulation for pout < 1e-5, as computation time would be very high
                if pout > 1e-5
                    NsamplesAdjusted = max( round( factor_sim * ( 1 / pout ) ), min_samples );
                    pout_sim_blocks(kg,ku,km,kfama) = SimOutage_BlocksFAMA(NsamplesAdjusted, gam_sim(kg), U(ku), sqrt(deltab), L, m(km), famatype(kfama), Inf);
                end
            end
        end
    end   
end
exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

figure(1)

% Colors for each "m" curve: "blue", "red" and "dark green"
colors = [[0 0 1]; [1 0 0]; [0.15,0.56,0.15]];

%---------------------------------------------------------------------
% Plot s-FAMA curves
%---------------------------------------------------------------------

subplot(2,1,1);
kfama = 1; 

for km = length(m):-1:1
    for ku = length(U):-1:1
         curveINT(km) = semilogy(gamdB, pout_calc_int_blocks(:,ku,km,kfama), 'Color', colors(km,:), 'linewidth', 3);
         hold on;
        curveGL = semilogy(gamdB_sim, pout_calc_gl_blocks(:,ku,km,kfama), 's', 'Color', 'r', 'linewidth', 3,'MarkerSize', 10);
        curveSIM = semilogy(gamdB_sim,  pout_sim_blocks(:,ku,km,kfama), '+', 'Color', 'k', 'linewidth', 3,'MarkerSize', 10);
    end
end


%-------------------------------------------------------------------------
% Set s-FAMA figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'SIR threshold $- \, \gamma$ (dB)', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [-15, 10] );
ylim( [0.99e-5, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curveINT(1), curveINT(2), curveINT(3), curveGL, curveSIM], ...
    {['$m=$ ',num2str(m(1))], ['$m=$ ',num2str(m(2))],['$m=$ ',num2str(m(3))],'Gauss-Laguerre', 'Simulation' },...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='southeast');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.479,0.612,0.21,0.075],...   
    'String',{'$s$-FAMA' , ...
    ['$N = $ ', num2str(N), ', $W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off');

% ------------ U = 10 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.288, 0.8, 0.0273, 0.0564]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.2167, 0.8285, 0.062, 0.0373],...
    'String','$U=10$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% ------------ U = 5 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.373, 0.7582, 0.0321, 0.069]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.41 0.76 0.062, 0.0373],...
    'String','$U=5$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%---------------------------------------------------------------------
% Plot f-FAMA curves
%---------------------------------------------------------------------

subplot(2,1,2);
kfama = 2; 

for km = length(m):-1:1
    for ku = length(U):-1:1
         curveINT(km) = semilogy(gamdB, pout_calc_int_blocks(:,ku,km,kfama), 'Color', colors(km,:), 'linewidth', 3);
         hold on;
        curveGL = semilogy(gamdB_sim, pout_calc_gl_blocks(:,ku,km,kfama), 's', 'Color', 'r', 'linewidth', 3,'MarkerSize', 10);
         hold on;
         curveSIM = semilogy(gamdB_sim,  pout_sim_blocks(:,ku,km,kfama), '+', 'Color', 'k', 'linewidth', 3,'MarkerSize', 10);
    end
end

%-------------------------------------------------------------------------
% Set f-FAMA figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'SIR threshold $- \, \gamma$ (dB)', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [-15, 10] );
ylim( [0.99e-5, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curveINT(1), curveINT(2), curveINT(3), curveGL, curveSIM], ...
    {['$m=$ ',num2str(m(1))], ['$m=$ ',num2str(m(2))],['$m=$ ',num2str(m(3))],'Gauss-Laguerre','Simulation' },...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='southeast');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.48, 0.147, 0.210, 0.075  ],...   
    'String',{'$f$-FAMA' , ...
    ['$N = $ ', num2str(N), ', $W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off');

% ------------ U = 10 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.288,0.268,0.0273,0.0799]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.267,0.353,0.062,0.0373],...
    'String','$U=10$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% ------------ U = 5 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.4598,0.3067,0.02969,0.069]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.503,0.292,0.062,0.0373],...
    'String','$U=5$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

