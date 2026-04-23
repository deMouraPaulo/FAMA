%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves for Figure 3.5
% Follow the steps:
% 1) Run the code Fig3_5.m to generate data; 
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



