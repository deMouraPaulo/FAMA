%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the Multiplexing Gain of fast O-FAMA network in terms 
% of the SIR threshold for a 1D fluid antenna under different correlation models: 
% 1) Simulation Jakes's (real), 
% 2) Proposed block diagonal expression (Eq. (29)), 
% 3) Simulation blocks, 
% 4) Average constant model in [16], 
% 5) Multiplexing gain approximation = min (U, M (1- p))
% MAKES FIGURE 2 IN THE PAPER
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
famatype = "Fast";          % Fast FAMA
N = 500;                     % Number of ports

W = 2;                       % Antenna size (wavelength normalized)
U = 50;                   % Number of users
M = 70;                     % Number of users in the pool for selection

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

gamdB = linspace(-10,20,20);     % SIR threshold (dB)
gam = 10.^(gamdB/10);            % SIR threshold (linear scale)

m = [1,2,3];                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------
% Jakes's model
mg_jakes = zeros(length(U),length(gam));


% Block-diagonal approximation (Lemma 2)
mg_blocks = zeros(length(U),length(gam));

% Approximated block-diagonal (Corollary 2)
mg_blocks_sim = zeros(length(U),length(gam));

% i.i.d. system
mg_approx = zeros(length(U),length(gam));

% Average correlation model in [16]
mg_avg = zeros(length(U),length(gam));

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


% Loop over fading (m)
for ku = 1:length(m )

    % User feedback
    disp(['Iter ' num2str(ku) ' out of ' num2str(length(m))]);

    %-----------------------------------------------------------------
    % Multiplexing Gain calculation/simulation
    %----------------------------------------------------------------- 
    % Jakes's outage probability
    p = SimOutageFAMA(Nsamples/100, gam, Sigma_jakes, U,m(ku),famatype, Inf);
    mg_jakes(ku,:) = MultiplexGain(M,U,p);

    % Simulated blocks outage probability
     p = SimOutage_BlocksFAMA(Nsamples/100, gam, U, sqrt(mu2), L, m(ku),famatype, Inf);
     mg_blocks_sim(ku,:) = MultiplexGain(M,U,p);


    % Block outage probability -> Evaluated through quadrature expression
    % in Eq. (31)
    p = CalcOutageFAMA(gam, L, mu2, U, 'Integral', order, m(ku),famatype,'SIR', Inf);
      mg_blocks(ku,:) = MultiplexGain(M,U,p);

   % Multiplexing Gain approximation
     mg_approx(ku,:) = min (U, M * (1-p));



    % Average correlation model 
    p = CalcOutageFAMA(gam, N, mu_avg^2, U, 'Integral',order, m(ku),famatype, 'SIR', Inf);
     mg_avg(ku,:) = MultiplexGain(M,U,p);


  
 

end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for ku = length(m):-1:1

    plot(gamdB, mg_jakes(ku,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    plot(gamdB, mg_blocks(ku,:), 'b-o', 'linewidth', 2);
    plot(gamdB, mg_approx(ku,:), ':','Marker','square', 'Color', [0.9290 0.6940 0.1250],'linewidth', 2);
    plot(gamdB, mg_blocks_sim(ku,:), 'r--^', 'linewidth', 2);
    plot(gamdB, mg_avg(ku,:), '--v', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2);
end

l = legend("Sim. Jakes", "Blocks", "Approximation" , "Sim. Blocks", "Constant");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('SIR threshold $\gamma$ (dB)', 'FontSize', 20, 'Interpreter','latex');
ylabel('Multiplexing Gain', 'FontSize', 20, 'Interpreter','latex');
ylim([0, U])
str = [num2str(famatype),' O-FAMA parameters: $N= $ ',num2str(N),', $W =$ ',num2str(W),', $U = $ ',num2str(U),', $M = $ ',num2str(M)];
title(str,'FontSize', 20, 'Interpreter','latex')



