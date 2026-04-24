%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code DON'T compute data base, ONLY plots curves for Figure 4.1
% Follow the steps:
% 1) Run the code Fig4_1.m to generate data; 
% 2) Run this code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)


colors = ['b';'r'];

kw = 1;
for kfama = length(famatype):-1:1
    curveINT(kfama) = semilogy(N, squeeze(pout_calc_int(:,kw,kfama)), 'Color', colors(kfama,:), 'linewidth', 4);
    hold on;
    curveINT_blocks(kfama) = semilogy(N, squeeze(pout_calc_int_blocks(:,kw,kfama)), 'Color', colors(kfama,:), 'linewidth', 4);
    curveGL = semilogy(Nsimul,squeeze(pout_calc_gl(:,kw,kfama)), 's', 'Color', 	"#77AC30", 'linewidth', 4,'MarkerSize', 10);
    curveGL_blocks = semilogy(Nsimul,squeeze(pout_calc_gl_blocks(:,kw,kfama)), 's', 'Color', 	"#77AC30", 'linewidth', 4,'MarkerSize', 10);
    curveSIM = semilogy(Nsimul, squeeze(pout_sim(:,kw,kfama)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    curveSIM_blocks = semilogy(Nsimul, squeeze(pout_sim_blocks(:,kw,kfama)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    curveSIM_jakes(kfama) = semilogy(Nsimul, squeeze(pout_sim_jakes(:,kw,kfama)), '--', 'Color', colors(kfama,:), 'linewidth', 2,'MarkerSize', 10);     
end



%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Ports $- \, N$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [N(1), N(end)] );
ylim( [6e-6, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT_blocks(1),curveSIM_jakes(1),curveINT_blocks(2),curveSIM_jakes(2), curveSIM, curveGL],...
    { '$s$-FAMA - calc. integ.', '$s$-FAMA - sim. Jakes','$f$-FAMA - calc. integ.','$f$-FAMA - sim. Jakes','Simulation SIR-based','Calc. Gauss-Laguerre'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='east');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.604363843034047,0.128860489882854,0.168328166275496,0.162517571884986],...   
    'String',['$W=$ ',num2str(W(kw)),newline, ...
    '$U=$ ',num2str(U), newline ...
    '$m = $ ', num2str(m), newline,...
    '$\gamma= $ ',num2str(gamdB),' dB'],...
    'Interpreter','latex',...
    'BackgroundColor','w',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','off');


%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Spatial-block correlation:
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.486583397982932,0.657108626198089,0.0394,0.252]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.361366408068269,0.758138019169338,0.125,0.074],...
    'String',{'Block Correlation'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------



%------------------------------------------------------------------------
% Const correlation:
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.149595805795877,0.418530351437705,0.424,0.074835995740152]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.318034445306439,0.346112886048989,0.125,0.057912673056444],...
    'String',{'Constant Correlation'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1],...
    'EdgeColor',[1 1 1]);
%------------------------------------------------------------------------

% ----------------------- FIGURE 2 --------------------------------------

figure(2)
colors = ['b';'r'];
kw = 2;
for kfama = length(famatype):-1:1
    curveINT(kfama) = semilogy(N, squeeze(pout_calc_int(:,kw,kfama)), 'Color', colors(kfama,:), 'linewidth', 4);
    hold on;
    curveINT_blocks(kfama) = semilogy(N, squeeze(pout_calc_int_blocks(:,kw,kfama)), 'Color', colors(kfama,:), 'linewidth', 4);
    curveGL = semilogy(Nsimul,squeeze(pout_calc_gl(:,kw,kfama)), 's', 'Color', 	"#77AC30", 'linewidth', 4,'MarkerSize', 10);
    curveGL_blocks = semilogy(Nsimul,squeeze(pout_calc_gl_blocks(:,kw,kfama)), 's', 'Color', 	"#77AC30", 'linewidth', 4,'MarkerSize', 10);
    curveSIM = semilogy(Nsimul, squeeze(pout_sim(:,kw,kfama)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    curveSIM_blocks = semilogy(Nsimul, squeeze(pout_sim_blocks(:,kw,kfama)), '+', 'Color', 'k', 'linewidth', 4,'MarkerSize', 10);
    curveSIM_jakes(kfama) = semilogy(Nsimul, squeeze(pout_sim_jakes(:,kw,kfama)), '--', 'Color', colors(kfama,:), 'linewidth', 2,'MarkerSize', 10);     
end



%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Ports $- \, N$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
xlim( [N(1), N(end)] );
ylim( [6e-6, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveINT_blocks(1),curveSIM_jakes(1),curveINT_blocks(2),curveSIM_jakes(2), curveSIM, curveGL],...
    { '$s$-FAMA - calc. integral', '$s$-FAMA - sim. Jakes','$f$-FAMA - calc. integral','$f$-FAMA - sim. Jakes','Simulation SIR-based','Calc. Gauss-Laguerre'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='southeast');

% Create textbox: parameters
annotation(figure(2),'textbox',...
    [0.666427458239632,0.517571884984028,0.168328166275496,0.162517571884987],...   
    'String',['$W=$ ',num2str(W(kw)),newline, ...
    '$U=$ ',num2str(U), newline ...
    '$m = $ ', num2str(m), newline,...
    '$\gamma= $ ',num2str(gamdB),' dB'],...
    'Interpreter','latex',...
    'BackgroundColor','w',...
    'HorizontalAlignment','center',...
    'FontSize',24,...
    'FitBoxToText','off');


%-------------------------------------------------------------------------
% Annotion boxes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Spatial-block correlation:
% Create ellipse
annotation(figure(2),'ellipse',...
    [0.489686578743211,0.405750798722045,0.0394,0.433439829605968]);
% Create textbox
annotation(figure(2),'textbox',...
    [0.359039022498059,0.654836634717796,0.125,0.074],...
    'String',{'Block Correlation'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%------------------------------------------------------------------------



%------------------------------------------------------------------------
% Const correlation:
% Create ellipse
annotation(figure(2),'ellipse',...
    [0.153631318305965,0.210862619808313,0.310452620354244,0.074835995740152]);
% Create textbox
annotation(figure(2),'textbox',...
    [0.197776172505793,0.322683706070289,0.125,0.057912673056444],...
    'String',{'Constant Correlation'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1],...
    'EdgeColor',[1 1 1]);
%------------------------------------------------------------------------


