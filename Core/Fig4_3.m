%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the Multipexig gain vs Number of users.
% Fast, Slow and Opportunistic FAMA networks with block-correlation model.
% 
% Makes Figure 4.3
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

% Number of users - slow FAMA ; Number of users - fast FAMA
U = round ([linspace(2,15,14); linspace(10,150,14)]); % (famatype, U) 

% Pool of users (famatype, U,  M)
Mpool = zeros(size(U,1),size(U,2),3);
Mpool(:,:,1) = round (1.5 * U);
Mpool(:,:,2) = round (2.0 * U);
Mpool(:,:,3) = round (2.5 * U);

% Transpose U (famatype, U )
% U = U';


    
gamdB = -3;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)


m = 2;                    % Nakagami-m fading severity
order = 50;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

mult_gain = zeros( size(U,1), size(U,2) ); 
mult_gain_OFAMA = zeros( size(U,1), size(U,2), size(Mpool,3) );


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

num_iter = numel(Mpool);
kiter = 0;


%-----------------------------------------------------------------
% Loop over  famatype
%-----------------------------------------------------------------
for kfama = 1:length(famatype)

    %-----------------------------------------------------------------
    % Loop over number of Users
    %-----------------------------------------------------------------
    for ku = 1: size (U,2) 
     
        % OP - calculation
        pout = CalcOutageFAMA(gam, L, deltab, U (kfama, ku) , 'Quadrature', order, m, famatype(kfama), 'SIR', Inf);

        % Multiplexing gain: s-FAMA and f-FAMA
        mult_gain(kfama, ku) = U(kfama, ku) * (1 - pout);

        %-----------------------------------------------------------------
        % Multiplexing gain O-FAMA
        % Loop over  Mpool
        %-----------------------------------------------------------------
        for kpool = 1:size(Mpool,3)

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            mult_gain_OFAMA(kfama, ku, kpool) = MultiplexGain(Mpool(kfama, ku, kpool), U(kfama, ku), pout );

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

% Creates a tiled chart layout for displaying multiple plots in the current figure
t = tiledlayout(2,1);

% ylabel for the two plots
ylabel(t, 'Multiplexing Gain $- \, \mathcal{G}_m $', 'Interpreter', 'Latex', 'FontSize', 24 );

% Tile 1
nexttile;

kfama = 1; 
for ku = size (U,2):-1:1
    curve = plot( U(kfama,:) , mult_gain(kfama, :),'-.', 'Color' ,'k', 'linewidth', 3 );
    hold on;
    for kpool = size(Mpool,3):-1:1
            curveOFAMA(kpool) = plot( U(kfama,:) , mult_gain_OFAMA(kfama, :,kpool), ...
                's-', 'Color', colors(kpool,:), 'linewidth', 3, 'MarkerSize', 10 );
    end
end


%-------------------------------------------------------------------------
% Set s-FAMA figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Users $- \, U$', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [2, 15] );
ylim( [0, 6] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curve, curveOFAMA(1), curveOFAMA(2), curveOFAMA(3)], ...
        {'$s$-FAMA','O-FAMA $(M = 1.5 \times U)$','O-FAMA $(M = 2 \times U)$','O-FAMA $(M = 2.5 \times U)$'},...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='northeast');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.1405,0.596,0.21,0.075],...   
    'String',{['$N = $, ', num2str(N), '$m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off');


%---------------------------------------------------------------------
% Plot f-FAMA curves
%---------------------------------------------------------------------

% Tile 2
nexttile;

kfama = 2; 
for ku = size (U,2):-1:1
    curve = plot( U(kfama,:) , mult_gain(kfama, :),'-.', 'Color' ,'k', 'linewidth', 3 );
    hold on;
    for kpool = size(Mpool,3):-1:1
            curveOFAMA(kpool) = plot( U(kfama,:) , mult_gain_OFAMA(kfama, :,kpool), ...
                's-', 'Color', colors(kpool,:), 'linewidth', 3, 'MarkerSize', 10 );
    end
end


%-------------------------------------------------------------------------
% Set f-FAMA figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Users $- \, U$', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [10, 150] );
ylim( [0, 150] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curve, curveOFAMA(1), curveOFAMA(2), curveOFAMA(3)], ...
        {'$s$-FAMA','O-FAMA $(M = 1.5 \times U)$','O-FAMA $(M = 2 \times U)$','O-FAMA $(M = 2.5 \times U)$'},...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='northwest');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.4256,0.358,0.21,0.075],...   
    'String',{['$N = $, ', num2str(N), '$m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off');

