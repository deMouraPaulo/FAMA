%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves for Figure 4.2
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


