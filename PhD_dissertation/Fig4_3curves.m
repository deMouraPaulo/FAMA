%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ONLY plot curves for Figure 4.3 (two figures in different files)
% Follow the steps:
% 1) Run the code Fig4_3.m to generate data; 
% 2) Run this code.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 



% Colors for each "m" curve: "blue", "red" and "dark green"
colors = [[0 0 1]; [1 0 0]; [0.15,0.56,0.15]];

%---------------------------------------------------------------------
% Plot s-FAMA curves
%---------------------------------------------------------------------

figure(1)

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
ylabel('Multiplexing Gain $- \, \mathcal{G}_m $', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [2, 15] );
ylim( [0, 6] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curve, curveOFAMA(1), curveOFAMA(2), curveOFAMA(3)], ...
        {'$s$-FAMA','O-FAMA $(M = 1.5 \times U)$','O-FAMA $(M = 2 \times U)$','O-FAMA $(M = 2.5 \times U)$'},...
     'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='southwest');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.6488,0.7763,0.2302,0.1193],'String',...   
    {'$s$-FAMA',...
    ['$N = $ ', num2str(N), '$, m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)]},...
    'BackgroundColor','w',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','off');


%---------------------------------------------------------------------
% Plot f-FAMA curves
%---------------------------------------------------------------------

figure(2)

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
ylabel('Multiplexing Gain $- \, \mathcal{G}_m $', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [10, 150] );
ylim( [0, 150] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curve, curveOFAMA(1), curveOFAMA(2), curveOFAMA(3)], ...
        {'$f$-FAMA','O-FAMA $(M = 1.5 \times U)$','O-FAMA $(M = 2 \times U)$','O-FAMA $(M = 2.5 \times U)$'},...
     'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='northeast');


% Create textbox: parameters
annotation(figure(2),'textbox',...
    [0.158293670886076,0.77310511182109,0.2302,0.11933],...   
    'String',{'$f$-FAMA',['$N = $ ', num2str(N), ', $m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'BackgroundColor','w',...
    'FontSize',24,...
    'FitBoxToText','off');


