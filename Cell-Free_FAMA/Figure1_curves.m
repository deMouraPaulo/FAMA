%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% This code DON'T compute data base, ONLY plots curves of a figure
% Follow the steps:
% 1) Run the code Figure1.m to generate data; 
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

for ka = 1:length(antBS)
semilogy(gamdB, pout_calcslow_int(ka,:), '-v', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2);
hold on; grid on;
semilogy(gamdB, pout_calcslow_quad(ka,:), '*', 'Color', 'b', 'linewidth', 2, MarkerSize=10);
semilogy(gamdB, pout_simCFslow(ka,:), '+', 'Color', 'k', 'linewidth', 2, MarkerSize=10);
semilogy(gamdB, pout_calcfast_int(ka,:), 'r.-','Marker','square', 'linewidth', 2);
semilogy(gamdB, pout_calcfast_quad(ka,:), '*', 'Color', 'b', 'linewidth', 2, MarkerSize=10);
semilogy(gamdB, pout_simCFfast(ka,:), '+','Color', 'k', 'linewidth', 2, MarkerSize=10)
end

leg = legend( "$s$-FAMA integral", "$s$-FAMA quadrature",  "$s$-FAMA simulation", ...
    "$f$-FAMA integral", "$f$-FAMA quadrature", "$f$-FAMA simulation",Location="northwest");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 
set(leg, 'FontSize', 24, 'Interpreter','latex');
xlabel('SIR threshold $\gamma$ (dB)', 'FontSize', 24, 'Interpreter','latex');
ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 24, 'Interpreter','latex');
xlim ([5, 20]);
ylim([1e-6, 1])


% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.330696202531644,0.14164004259851,0.233386075949368,0.106766847693775],...   
    'String',{['$N = $ ', num2str(N), ', $W = $ ',num2str(W)], ...
    ['$U = $ ', num2str(U), ',  m = ',num2str(m)], ...
    ['$d_0=$ ',num2str(d0), ', $d=$ ',num2str(d), ', $\alpha = $ ',num2str(alphaCF) ] },...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'BackgroundColor','w',...
    'FontSize',20,...
    'FitBoxToText','off');


% ------------ nant = 2 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.547493670886075,0.608093716719915,0.0273,0.273695420660286]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.530782278481011,0.562259318423867,0.174913924050635,0.0373],...
    'String','BS antennas $=2$',...
    'Interpreter','latex',...
    'FontSize',20,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% ------------ nant = 8 ---------------------------
% Create ellipse
annotation(figure(1),'ellipse',...
    [0.719962025316454,0.42066027689031,0.0273,0.236421725239616]);
% Create textbox
annotation(figure(1),'textbox',...
    [0.709579746835438,0.365241214057523,0.174122784810131,0.0373],...
    'String','BS antennas $=8$',...
    'Interpreter','latex',...
    'FontSize',20,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'EdgeColor','none');



%---------------------- Mean and Variance Analysis ------------------------

%-------- Nakagmi-m  RVs --------------------------------------------------

% sinal desejado
m_des = m .* antBS;
omega_des =  2 .* m .* antBS .* d0.^(-alphaCF);

mean_des = ( (gamma(m_des + 0.5) ./ gamma(m_des)) ) .* sqrt(omega_des ./ m_des) ;
var_des = omega_des .* ( 1 - (1./m_des) .* ( gamma(m_des + 0.5) ./ gamma(m_des) ).^2 );

% interf s-FAMA
m_sFAMA = m .* U ;
omega_sFAMA = 2 .* m .* U .* d.^(-alphaCF);
%
mean_slow = ( (gamma(m_sFAMA + 0.5) ./ gamma(m_sFAMA)) ) .* sqrt(omega_sFAMA ./ m_sFAMA) ;
var_slow = omega_sFAMA .* ( 1 - (1./m_sFAMA) .* ( gamma(m_sFAMA + 0.5) ./ gamma(m_sFAMA) ).^2 );

% interf f-FAMA
m_fFAMA = (m .* U) ./ (1 + m.*(U-1));
omega_fFAMA = omega_sFAMA;
% 
mean_fast = ( (gamma(m_fFAMA + 0.5) ./ gamma(m_fFAMA)) ) .* sqrt(omega_fFAMA ./ m_fFAMA) ;
var_fast = omega_fFAMA .* ( 1 - (1./m_fFAMA) .* ( gamma(m_fFAMA + 0.5) ./ gamma(m_fFAMA) ).^2 );

disp('NAKAGAMI-M RVs:')
disp(['Mean of the desired envelope: ', num2str(mean_des) ])
disp(['Variance of the desired envelope: ', num2str(var_des) ])
disp('---')
disp(['Ratio of the means desired: ', num2str(mean_des(2)/mean_des(1)) ])
disp(['Ratio of the variances desired: ', num2str(var_des(2)/var_des(1)) ])
disp('---')
disp(['Mean of the interference envelope, f-FAMA: ', num2str(mean_fast) ])
disp(['Mean of the interference envelope, s-FAMA: ', num2str(mean_slow) ])
disp(['Variance of the interference envelope, f-FAMA: ', num2str(var_fast) ])
disp(['Variance of the interference envelope, s-FAMA: ', num2str(var_slow) ])
disp('----')
disp(['Ratio of the means (interference envelope), f-FAMA/s-FAMA: ', num2str(mean_fast/mean_slow) ])
disp(['Ratio of the variances (interference envelope), f-FAMA/s-FAMA: ', num2str(var_fast/var_slow) ])

%------------------- chi-square (X^2) RVs --------------------------------------

% mean ( X^2 ) = n * sig^2
% where: n: number of degrees of freedom (DoFs), and 
% sig^2 = variance of the Gaussian RV that are added to form the X^2 RV
%
% variance ( X^2 ) = 2 * n * ( sig^2 )^2
%
% For desired user
%
% DoFs = 2 * m * antBS
% var gauss =  d0.^(-alphaCF)
%
mean_desX2 = 2 .* m .* antBS .* d0.^(-alphaCF);
var_desX2 = 4 .* m .* antBS .* d0.^(-2*alphaCF);
%
% interf s-FAMA
%
% DoFs = 2 * m * U
% var gauss =  d0.^(-alphaCF)
%
mean_slowX2 = 2 .* m .* U .* d.^(-alphaCF);
var_slowX2 = 4 .* m .* U .* d.^(-2*alphaCF);
%
% interf f-FAMA
%
% DoF = 2
% Var Gaussian = d.^(-alphaCF) * Util / mtil = d.^(-alphaCF) m * U * (m*(U-1)+1) / m * U) 
% = ( m * (U-1) + 1) * d.^(-alphaCF)
% obs.: para mtil = 1 -> Var Gauss = Util = m * U * d.^(-alphaCF)
%
mean_fastX2 = 2 .* (m .* (U-1) + 1) .* d.^(-alphaCF);
var_fastX2 = 4 .* ( (m .* (U-1) + 1) .* d.^(-alphaCF) )^2;


disp('---')
disp('---')
disp('Chi-square RVs')
disp(['Mean of the desired envelope: ', num2str(mean_desX2) ])
disp(['Variance of the desired envelope: ', num2str(var_desX2) ])
disp('---')
disp(['Ratio of the means desired: ', num2str(mean_desX2(2)/mean_desX2(1)) ])
disp(['Ratio of the variances desired: ', num2str(var_desX2(2)/var_desX2(1)) ])
disp('---')
disp(['Mean of the interference envelope, f-FAMA: ', num2str(mean_fastX2) ])
disp(['Mean of the interference envelope, s-FAMA: ', num2str(mean_slowX2) ])
disp(['Variance of the interference envelope, f-FAMA: ', num2str(var_fastX2) ])
disp(['Variance of the interference envelope, s-FAMA: ', num2str(var_slowX2) ])
disp('----')
disp(['Ratio of the means (interference envelope), f-FAMA/s-FAMA: ', num2str(mean_fastX2/mean_slowX2) ])
disp(['Ratio of the variances (interference envelope), f-FAMA/s-FAMA: ', num2str(var_fastX2/var_slowX2) ])
