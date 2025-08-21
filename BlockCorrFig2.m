%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP in terms of the SIR threshold for a 1D fluid
% antenna under different correlation models: 
% 1) Simulation Jakes's (real), 
% 2) Proposed block diagonal expression (Eq. (29)), 
% 3) Simulation blocks, 
% 4) Average constant model, 
% 5) Upper bound = i.i.d. (Eq. (38))
% MAKES FIGURE 1 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------
clc
close all
clear

addpath('Core/')

%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------

famatype = "Slow";

N = 100;                     % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = 5;                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

gamdB = linspace(-10,10,20);     % SIR threshold (dB)
gam = 10.^(gamdB/10);            % SIR threshold (linear scale)

m = [1,2,3];                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------
% Jakes's model
pout_jakes = zeros(length(U),length(gam));

% Block-diagonal approximation (Lemma 2)
pout_blocks = zeros(length(U),length(gam));

% Approximated block-diagonal (Corollary 2)
pout_blocks_sim = zeros(length(U),length(gam));

% i.i.d. system
pout_iid = zeros(length(U),length(gam));

% Average correlation model in [16]
pout_avg = zeros(length(U),length(gam));

%-----------------------------------------------------------------
% Jake's correlation
%----------------------------------------------------------------
Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));

rho = sort(eig(Sigma_jakes),'descend');

%----------------------------------------------------------------
% Block diagonal correlation matrix approximation
%-----------------------------------------------------------------
mu2 = 0.97;
Num_eig = sum(rho > N/100);

% Algorithm 1
L = BlockCorrelation(N, rho, Num_eig, mu2);

% Average correlation coefficient in [16]
mu_avg = sqrt(2)*sqrt(hypergeom(0.5, [1 1.5], -pi^2*W^2) - ...
                besselj(1, 2*pi*W)/(2*pi*W));


% Loop over number of users (m)
for ku = 1:length(m )

    % User feedback
    disp(['Iter ' num2str(ku) ' out of ' num2str(length(m))]);

    %-----------------------------------------------------------------
    % Outage Probabilities calculation/simulation
    %----------------------------------------------------------------- 
    % Jakes's outage probability
    pout_jakes(ku,:) = SimOutageFAMA(Nsamples, gam, Sigma_jakes, U,m(ku),famatype,Inf);

    % Simulated blocks outage probability
    pout_blocks_sim(ku,:) = SimOutage_BlocksFAMA(Nsamples,gam, U, sqrt(mu2), L, m(ku),famatype, Inf);

    % Block outage probability -> Evaluated through quadrature expression in Eq. (31) 
     pout_blocks(ku,:) = CalcOutageFAMA(gam, L, mu2, U, 'Quadrature', order, m(ku),famatype,'SIR',Inf);

    % Average correlation model in [16]  
    pout_avg(ku,:) = CalcOutageFAMA(gam, N, mu_avg^2, U, 'Quadrature',order, m(ku),famatype,'SIR',Inf);


    % I.I.D. outage =  Upper bound outage
%     pout_iid(ku,:) = (1 - 1./(1+gam).^(U(ku)-1)).^Num_eig; 
    pout_iid(ku,:) = OutageUBblocks(gam, Num_eig, U, m(ku),famatype);

  

end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for ku = length(m):-1:1

    semilogy(gamdB, pout_jakes(ku,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(gamdB, pout_blocks(ku,:), 'b-o', 'linewidth', 2);
    semilogy(gamdB, pout_avg(ku,:), '-v', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2);
    semilogy(gamdB, pout_iid(ku,:), 'r.-','Marker','square', 'linewidth', 2);
    semilogy(gamdB, pout_blocks_sim(ku,:), 'c--^', 'Color', [0.9290 0.6940 0.1250], 'linewidth', 2);
end

l = legend( "Simulation", "Blocks", "Constant", "Upper Bound" ,"Blocks Sim.");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('SIR threshold $\gamma$ (dB)', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{\mathrm{out}}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-5, 1])
str = [num2str(famatype),' FAMA parameters: $N= $ ',num2str(N),', $W =$ ',num2str(W),', $U = $ ',num2str(U)];
title(str,'FontSize', 20, 'Interpreter','latex')



