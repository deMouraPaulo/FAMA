%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves for Figure 3.4
% Follow the steps:
% 1) Run the code Fig3_4.m to generate data; 
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

% Colors for each "m" curve: "light blue", "light organge and "red"";
colors = [ [71,147,175]; [ 255, 196, 110];  [255, 0, 0] ]/255;

for km = length(m):-1:1
    for ku = length(U):-1:1
        curveINT(km) = loglog(N, squeeze(pout_calc_int(:,km,ku)), 'Color', colors(km,:), 'linewidth', 4);
        hold on;
        curveGL = loglog(Nsimul,squeeze(pout_calc_gl(:,km,ku)), 's', 'Color', 'r', 'linewidth', 4,'MarkerSize', 10);
        curveSIM = loglog(Nsimul, squeeze(pout_sim(:,km,ku)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    end
end

%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Ports $- \, N$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [10, 150] );
ylim( [0.99e-5, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT(1), curveINT(2), curveINT(3), curveGL, curveSIM],...
    { ['$m = $ ', num2str(m(1))], ['$m= $ ',num2str(m(2))], ['$m= $ ',num2str(m(3))],'Gauss-Laguerre','Simulation'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='west');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.464,0.133,0.3161,0.1471],...   
    'String',['Average SNR $=$ ',num2str(avg_snrdB), ' dB', newline ...
    '$\gamma = $ ', num2str(gamdB), ' dB' newline,...
    '$W =$ ',num2str(W)],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'BackgroundColor','w',...
    'FontSize',24,...
    'FitBoxToText','off');


%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% U = 5
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.7983,0.720,0.0394,0.205]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.783,0.667,0.0826,0.0447],...
    'String',{'$U = 5$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% U = 3
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.540,0.460,0.0394,0.316]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.511,0.405,0.0826,0.0447],...
    'String',{'$U = 3$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------



