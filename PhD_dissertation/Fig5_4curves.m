%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves of a figure
% Follow the steps:
% 1) Run the code Figure5_4.m to generate data; 
% 2) Run this code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

close all

%---------------------------------------------------------------------
% Plot multiplexing gain curves - figure 1
%---------------------------------------------------------------------

figure(1)

% plot s-fama
kfama = 1; 
curve_s = plot( U -1, mult_gain(kfama, :),'-s', 'Color' ,'r', 'linewidth', 3, 'MarkerSize', 10 );
hold on;
% plot f-fama
kfama = 2; 
curve_f = plot( U -1, mult_gain(kfama, :),'--s', 'Color' ,'b', 'linewidth', 3, 'MarkerSize', 10 );


%-------------------------------------------------------------------------
% Set multiplexing gain figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of interfering BS $- \, U$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Multiplexing Gain $- \, \mathcal{G}_m $', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [min(U-1), max(U)-1] );
ylim( [0, max(max(mult_gain))] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curve_s, curve_f ], ...
        {'$s$-FAMA', '$f$-FAMA'},...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='northwest');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.387565977798011,0.141640042598513,0.278012658227848,0.107634717784881],...   
    'String',{['$N = $ ', num2str(N),', K = ',num2str(antBS) ,'$, m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)], ...
    ['d0 = ',num2str(d0),', d = ',num2str(d), ', $\alpha = $ ',num2str(alphaCF)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'BackgroundColor','w',...
    'FontSize',20,...
    'FitBoxToText','off');


%---------------------------------------------------------------------
% Plot outage probability curves - figure 2
%---------------------------------------------------------------------

figure(2)

% plot s-fama
kfama = 1; 
curve_s = semilogy( U -1, pout(kfama, :),'-s', 'Color' ,'r', 'linewidth', 3, 'MarkerSize', 10 );
hold on;
% plot f-fama
kfama = 2; 
curve_f = semilogy( U -1, pout(kfama, :),'-s', 'Color' ,'b', 'linewidth', 3, 'MarkerSize', 10 );

%-------------------------------------------------------------------------
% Set OP figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of interfering BS $- \, U$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability ', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [min(U)-1, max(U)-1] );
ylim( [ min(min(pout)), 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curve_s, curve_f ], ...
        {'$s$-FAMA', '$f$-FAMA'},...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='northwest');


% Create textbox: parameters
annotation(figure(2),'textbox',...
    [0.387625316455696,0.130990415335464,0.257944303797468,0.110316293929713],...   
    'String',{['$N = $ ', num2str(N),', K = ',num2str(antBS) ,'$, m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)], ...
    ['d0 = ',num2str(d0),', d = ',num2str(d), ', $\alpha = $ ',num2str(alphaCF)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'BackgroundColor','w',...
    'FontSize',20,...
    'FitBoxToText','off');


