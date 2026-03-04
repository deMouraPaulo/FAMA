%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves of a figure
% Follow the steps:
% 1) Run the code Figure2.m to generate data; 
% 2) Run this code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

close all

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

figure(1)

% Colors for each "m" curve: "blue", "red" and "dark green"
colors = [[0 0 1]; [1 0 0]; [0.15,0.56,0.15]];

%---------------------------------------------------------------------
% Plot multiplexing gain curves
%---------------------------------------------------------------------

% Creates a tiled chart layout for displaying multiple plots in the current figure
t = tiledlayout(2,1);

% Tile 1
nexttile;

for ku = size (U,2):-1:1
    kfama = 1; 
    curve_s = plot( U -1, mult_gain(kfama, :),'--', 'Color' ,'r', 'linewidth', 3 );
    hold on;
    kfama = 2; 
    curve_f = plot( U -1, mult_gain(kfama, :),':', 'Color' ,'b', 'linewidth', 3 );
end


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
    [0.270246835443038,0.802981895633653,0.278012658227848,0.107634717784881],...   
    'String',{['$N = $ ', num2str(N),', K = ',num2str(antBS) ,'$, m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)], ...
    ['d0 = ',num2str(d0),', d = ',num2str(d), ', $\alpha = $ ',num2str(alphaCF)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'BackgroundColor','w',...
    'FontSize',20,...
    'FitBoxToText','off');


%---------------------------------------------------------------------
% Plot outage probability curves
%---------------------------------------------------------------------

% Tile 2
nexttile;

% for ku = size (U,2):-1:1
    kfama = 1; 
    curve_s = semilogy( U -1, pout(kfama, :),'-', 'Color' ,'r', 'linewidth', 3 );
    hold on;
    kfama = 2; 
    curve_f = semilogy( U -1, pout(kfama, :),'-', 'Color' ,'b', 'linewidth', 3 );
% end


%-------------------------------------------------------------------------
% Set OP figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of interfering BS $- \, U$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( '$P_{\mathrm{out}}(U)$ ', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [min(U)-1, max(U)-1] );
ylim( [ min(min(pout)), 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 

legend([curve_s, curve_f ], ...
        {'$s$-FAMA', '$f$-FAMA'},...
     'Interpreter', 'Latex', 'FontSize', 20, 'NumColumns',1, Location='northwest');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.387625316455696,0.130990415335464,0.257944303797468,0.110316293929713],...   
    'String',{['$N = $ ', num2str(N),', K = ',num2str(antBS) ,'$, m = $ ', num2str(m) ], ...
    ['$\gamma = $ ', num2str(gamdB), ' dB, ', '$W = $ ',num2str(W)], ...
    ['d0 = ',num2str(d0),', d = ',num2str(d), ', $\alpha = $ ',num2str(alphaCF)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'BackgroundColor','w',...
    'FontSize',20,...
    'FitBoxToText','off');


%------------------- chi-square (X^2) RVs --------------------------------------

% For desired user
%
mean_desX2 = 2 .* m .* antBS .* d0.^(-alphaCF);
var_desX2 = 4 .* m .* antBS .* d0.^(-2*alphaCF);

% interf s-FAMA
%
mean_slowX2 = 2 .* m .* U .* d.^(-alphaCF);
var_slowX2 = 4 .* m .* U .* d.^(-2*alphaCF);
%
% interf f-FAMA
%
mean_fastX2 = 2 .* (m .* (U-1) + 1) .* d.^(-alphaCF);
var_fastX2 = 4 .* ( (m .* (U-1) + 1) .* d.^(-alphaCF) ).^2;


disp('Chi-square RVs: ')
disp('Number of interfering BS:')
disp (num2str(U-1))
disp('Ratio of the means (interference envelope), f-FAMA/s-FAMA: ')
disp ( num2str(mean_fastX2./mean_slowX2))
disp('Ratio of the variances (interference envelope), f-FAMA/s-FAMA: ' )
disp(num2str(var_fastX2./var_slowX2))
disp('---')
disp(['Mean of the desired envelope: ', num2str(mean_desX2) ])
disp(['Variance of the desired envelope: ', num2str(var_desX2) ])
