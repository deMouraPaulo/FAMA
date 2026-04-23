%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots OP for FAS and MRC receivers
%
% Makes Figure 2.3
%
% For FAS OP, it was set U = 1 (one user) in the function CalcOutageFAMA.m
%
% For MRC OP, it used Eq. (6.25) from [1].
% [1] G. L. Stüber, Principles of Mobile Communication, 4th ed., 
% Cham, Switzerland: Springer International Publishing, 2017.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

tic
clc
clear
close all
addpath(fullfile('..','Core'))


%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------


N = [2,3,5,10,20,50,80,100,120,170];      % Number of ports
W = [0.2,0.5,1,2,5];                      % Antenna size (in normalized wavelength)

famatype = "Fast";    % FAMA type (for U = 1, don't care)
U = 1;                % Number of users (for FAS U = 1)
m = 1;                % Rayleigh

gam_avg = 1;
gam = 1;


% ant MRC
n_ant = [2,4,6,8,9];



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_FAS = zeros(length(W),length(N));
pout_calc_MRC = zeros(length(n_ant),1);


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

%------- Pout MRC -------------------
for kant = 1:(length(n_ant))
    soma = 0;
    for  k = 0:n_ant(kant)-1
       soma = soma + (1/factorial(k)) * (gam/gam_avg)^k;
    end
    pout_calc_MRC(kant) = 1 - exp(-gam/gam_avg)*soma;
end

%--------- Pout FAS ---------------------

% Loop over  W
for kw = 1:length(W)
    % User feedback
    disp(['Iter "W": ' num2str(kw) ' out of ' num2str(length(W))]); 
    % Loop over N
    for kn = 1:length(N)
        % Average correlation coefficient
        dm = DeltaMed (N(kn), W(kw));
        % Pout
        pout_calc_FAS(kw,kn) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral', 1, m,famatype,1,1, gam_avg);
    end
end


%---------------------------------------------------------------------
% Plotting
%---------------------------------------------------------------------

figure(1)

% Plot FAS
for kw = length(W):-1:1
    curveFAS(kw) = loglog(N, pout_calc_FAS(kw,:), '-*',  'linewidth', 4);
    hold on;
end

% Plot MRC
for kant = length(n_ant):-1:1
    curvMRC = loglog(N, pout_calc_MRC(kant)*ones(1,length(N)), '--k' ,'LineWidth',2);
end

grid on;
xlabel( 'Number of Ports $- \, N$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 24 );
 xlim( [min(N), max(N)] );
 ylim( [0.99e-6, 1] );
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveFAS(1), curveFAS(2), curveFAS(3), curveFAS(4), curveFAS(5), curvMRC],...
    { ['$W = $ ', num2str(W(1))], ['$W= $ ',num2str(W(2))], ['$W= $ ',num2str(W(3))], ...
    ['$W= $ ',num2str(W(4))],['$W= $ ',num2str(W(5))],'MRC'},...
    'Interpreter', 'Latex', 'FontSize', 24, 'NumColumns',1, Location='southwest');


% Num MRC Antennas = 9
% Create textbox
annotation(figure(1),'textbox',...
    [0.328688130333591,0.132369542066032,0.184888285492631,0.0447],...
    'String',{'Antennas = $9$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Num MRC Antennas = 8
% Create textbox
annotation(figure(1),'textbox',...
    [0.330239720713731,0.263359957401496,0.184888285492631,0.0447],...
    'String',{'Antennas = $8$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');


% Num MRC Antennas = 6
% Create textbox
annotation(figure(1),'textbox',...
    [0.252660201706748,0.500846645367419,0.184888285492631,0.0447],...
    'String',{'Antennas = $6$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');


% Num MRC Antennas = 4
% Create textbox
annotation(figure(1),'textbox',...
    [0.132411947245926,0.706384451544206,0.184888285492631,0.0447],...
    'String',{'Antennas = $4$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');


% Num MRC Antennas = 2
% Create textbox
annotation(figure(1),'textbox',...
    [0.690208688906127,0.773477103301399,0.184888285492631,0.0447],...
    'String',{'Antennas = $2$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','right',...
    'FontSize',24,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');
% Create line
annotation(figure(1),'line', ...
    [0.712179984484096 0.660977501939488],...
    [0.798787007454741 0.842385516506924], ...
    'LineWidth',1);

