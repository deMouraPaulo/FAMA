%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs Number of ports.
% Fast and slow FAMA networks, constant and block-correlation models.
% 
% Makes Figure 4.1.
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
N = 20:5:200;                   % Number of ports
Nsimul = [20,60,100,140,180];    % Number of ports for simulation
W = [1,3];                      % Antenna size (wavelength normalized)
U = 5;                          % Number of users

min_samples = 5e5;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e2;            % Multiplication factor for simulation 
     
gamdB = -3;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)


m = 2;                    % Nakagami-m fading severity
order = 50;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_int = zeros(length(N), length(W), length(famatype));
pout_calc_int_blocks = zeros(length(N), length(W), length(famatype));

pout_calc_gl = zeros(length(Nsimul), length(W), length(famatype));
pout_sim = zeros(length(Nsimul), length(W), length(famatype));

pout_calc_gl_blocks = zeros(length(Nsimul), length(W), length(famatype));
pout_sim_blocks = zeros(length(Nsimul), length(W), length(famatype));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter = length(N) *  length(W) * length(famatype);
kiter = 0;
% Loop over  W
for kw = 1:length(W)

    % User feedback
    disp(['Iter "W": ' num2str(kw) ' out of ' num2str(length(W))]);

    % Loop over  famatype
    for kfama = 1:length(famatype)

        %-----------------------------------------------------------------
        % Loop over N - Integral calculation
        %-----------------------------------------------------------------
        for kn = 1:length(N )

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            %-----------------------------------------------------------------
            % Constant Correlation
            %-----------------------------------------------------------------
            % Average correlation coefficient
            dm = DeltaMed (N(kn), W(kw));

            % Integral calculation
            pout_calc_int(kn,kw,kfama) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral', order, m, famatype(kfama), 'SIR', Inf);        
            
            %-----------------------------------------------------------------
            %  Block correlation
            %----------------------------------------------------------------
            % Correlation matriz
            Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N(kn)-1)*W(kw)/(N(kn)-1)));
            % Eigenvalues 
            rho = sort(eig(Sigma_jakes),'descend');
            % Correlation coefficent per block
            deltab = 0.97;
            % Number of domminant eigenvalues
            Num_eig = sum(rho > N(kn)/100);
            % Algorithm 1. L: vector with block sizes (Lb)
            L = BlockCorrelation(N(kn), rho, Num_eig, deltab);

            % Integral calculation
            pout_calc_int_blocks(kn,kw,kfama) = CalcOutageFAMA(gam, L, deltab, U, 'Integral', order, m, famatype(kfama), 'SIR', Inf);

        end

        %-----------------------------------------------------------------
        % Loop over Nsimul - Simulation and Quadrature Calculation
        %-----------------------------------------------------------------
        for kn = 1:length( Nsimul )

            %-----------------------------------------------------------------
            % Constant Correlation
            %-----------------------------------------------------------------
            % Average correlation coefficient
            dm = DeltaMed (Nsimul(kn), W(kw));

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl(kn,kw,kfama) = CalcOutageFAMA(gam, Nsimul(kn), dm, U, 'Quadrature', order, m,famatype(kfama), 'SIR', Inf);
            
            %------------ Monte Carlo simulation -------------------------

            % Pout calculated
            pout = pout_calc_gl(kn,kw,kfama);

            % Avoid simulation for pout < 0.9e-5, as computation time would be very high
            if pout > 0.9e-5
                NsamplesAdjusted = max( round( factor_sim * ( 1 / pout ) ), min_samples );
                pout_sim(kn,kw,kfama) = SimOutage_BlocksFAMA(NsamplesAdjusted, gam, U, sqrt(dm), Nsimul(kn), m,famatype(kfama), Inf);
            end
            %-----------------------------------------------------------------
            %  Block correlation
            %----------------------------------------------------------------
            % Correlation matriz
            Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:Nsimul(kn)-1)*W(kw)/(Nsimul(kn)-1)));
            % Eigenvalues 
            rho = sort(eig(Sigma_jakes),'descend');
            % Correlation coefficent per block
            deltab = 0.97;
            % Number of domminant eigenvalues
            Num_eig = sum(rho > Nsimul(kn)/100);
            % Algorithm 1. L: vector with block sizes (Lb)
            L = BlockCorrelation(Nsimul(kn), rho, Num_eig, deltab);

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl_blocks(kn,kw,kfama) = CalcOutageFAMA(gam, L, deltab, U, 'Quadrature', order, m,famatype(kfama), 'SIR', Inf);

            % Monte Carlo simulation
            NsamplesAdjusted = max( round( factor_sim * ( 1 / pout_calc_gl_blocks(kn,kw,kfama) ) ), min_samples );
            pout_sim_blocks(kn,kw,kfama) = SimOutage_BlocksFAMA(NsamplesAdjusted, gam, U, sqrt(deltab), L, m, famatype(kfama), Inf);
            
        end
    end   
end
exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

% Colors for SIR curves: "light orange" and "light blue";
colors = ['b';'r'];

for kw = length(W):-1:1
    for kfama = length(famatype):-1:1
        curveINT(kw) = semilogy(N, squeeze(pout_calc_int(:,kw,kfama)), 'Color', colors(kw,:), 'linewidth', 4);
        hold on;
        curveINT_blocks(kw) = semilogy(N, squeeze(pout_calc_int_blocks(:,kw,kfama)), 'Color', colors(kw,:), 'linewidth', 4);
        curveGL = semilogy(Nsimul,squeeze(pout_calc_gl(:,kw,kfama)), 's', 'Color', 'r', 'linewidth', 4,'MarkerSize', 10);
        curveGL_blocks = semilogy(Nsimul,squeeze(pout_calc_gl_blocks(:,kw,kfama)), 's', 'Color', 'r', 'linewidth', 4,'MarkerSize', 10);
        curveSIM = semilogy(Nsimul, squeeze(pout_sim(:,kw,kfama)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
        curveSIM_blocks = semilogy(Nsimul, squeeze(pout_sim_blocks(:,kw,kfama)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    end
end

%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Ports $- \, N$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [20, 200] );
ylim( [6e-6, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT(1), curveINT(2), curveGL, curveSIM],...
    { ['$W = $ ', num2str(W(1))], ['$W = $ ',num2str(W(2))], 'Gauss-Laguerre','Simulation'},...
    'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='east');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.627, 0.128, 0.112, 0.108],...   
    'String',['$U=$ ',num2str(U), newline ...
    '$m = $ ', num2str(m), newline,...
    '$\gamma= $ ',num2str(gamdB),' dB'],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off');


%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Spatial-block correlation s-FAMA:
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.4558, 0.760, 0.0394, 0.133]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.4917, 0.8444, 0.125, 0.074],...
    'String',{'Block Corr.','$s$-FAMA'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Spatial-block correlation f-FAMA:
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.454, 0.475, 0.0394, 0.252]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.355, 0.566, 0.125, 0.074],...
    'String',{'Block Corr.','$f$-FAMA'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Const correl s-FAMA:
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.263,0.34,0.424,0.049]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.3576,0.346,0.125,0.074],...
    'String',{'Cte. Corr.','$s$-FAMA'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1],...
    'EdgeColor',[1 1 1]);
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Const correl f-FAMA:
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.1306,0.1629,0.221,0.113]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.1855,0.175,0.124,0.074],...
    'String',{'Cte. Corr.','$f$-FAMA'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);
%------------------------------------------------------------------------


