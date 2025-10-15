%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves for Figure 4.3
% Follow the steps:
% 1) Run the code Fig4_3.m to generate data; 
% 2) Run this code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close

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


