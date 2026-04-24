%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves for Figure 4.2
% Figure 4.2 (two figures in different files)
% Follow the steps:
% 1) Run the code Fig4_2.m to generate data; 
% 2) Run this code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

figure(1)

% Colors for each "m" curve: "blue", "red" and "dark green"
colors = [[0 0 1]; [1 0 0]; [0.15,0.56,0.15]];

%---------------------------------------------------------------------
% Plot s-FAMA curves
%---------------------------------------------------------------------

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
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT(1), curveINT(2), curveINT(3), curveGL, curveSIM], ...
    {['$m=$ ',num2str(m(1))], ['$m=$ ',num2str(m(2))],['$m=$ ',num2str(m(3))],'Gauss-Laguerre', 'Simulation' },...
     'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='southeast');


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.442025316455696,0.18853354632588,0.21,0.075],...   
    'String',{'$s$-FAMA' , ...
    ['$N = $ ', num2str(N), ', $W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','off');

% ------------ U = 10 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.373443037974683,0.783812566560174,0.0273,0.119445793397234]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.273662025316456,0.847669329073487,0.062,0.0373],...
    'String','$U=10$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% ------------ U = 5 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.45686075949367,0.743343982960599,0.0321,0.099830457933977]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.506256962025316,0.761407348242816,0.062,0.0373],...
    'String','$U=5$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%---------------------------------------------------------------------
% Plot f-FAMA curves
%---------------------------------------------------------------------

figure(2)
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
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT(1), curveINT(2), curveINT(3), curveGL, curveSIM], ...
    {['$m=$ ',num2str(m(1))], ['$m=$ ',num2str(m(2))],['$m=$ ',num2str(m(3))],'Gauss-Laguerre','Simulation' },...
     'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='southeast');

% Create textbox: parameters
annotation(figure(2),'textbox',...
    [0.442025316455696,0.18853354632588,0.21,0.075 ],...   
    'String',{'$f$-FAMA' , ...
    ['$N = $ ', num2str(N), ', $W = $ ',num2str(W)]},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','off');

% ------------ U = 10 ---------------------------
% Create ellipse
annotation(figure(2),'ellipse',...
    [0.288,0.490947816826412,0.0273,0.141297231096915]);
% Create textbox
annotation(figure(2),'textbox',...
    [0.200544303797468,0.609656017039406,0.062,0.0373],...
    'String','$U=10$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% ------------ U = 5 ---------------------------
% Create ellipse
annotation(figure(2),'ellipse',...
    [0.458217721518987,0.582534611288609,0.02969,0.130758572949951]);
% Create textbox
annotation(figure(2),'textbox',...
    [0.502208860759494,0.602969116080938,0.062,0.0373],...
    'String','$U=5$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');


