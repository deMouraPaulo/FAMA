%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs SIR threshold gamma(dB)
% Fast and slow FAMA networks with block-correlation model.
% 
% Makes Figure 4.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

tic
clc
clear
close all
addpath(fullfile('..','Core'))



%-------------------------------------------------------------------------
% Flags
%-------------------------------------------------------------------------

% set flags for f-FAMA calculus (for s-Fama don't care)

% flag, approx. mf = 1
flagmf = true;

% flag, 2nd integral
flag2int = true;

% flag, simulation f-FAMA, 
flagsimul = true; % (mf = 1, approx. method, Gaussian RVs)


%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------

famatype = ["Slow";"Fast"];     %  FAMA type
N = 100;                   % Number of ports
W = 1;                      % Antenna size (wavelength normalized)
U = [5, 10];                % Number of users

min_samples = 1e4;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e2;            % Multiplication factor for simulation 
     
gamdB = -15:10;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)

gamdB_sim =  linspace(-15,10,10);   % SIR threshold (dB) for simulation
gam_sim = db2pow(gamdB_sim);        % SIR threshold (linear scale) for simulation


m = [1,2,3];                    % Nakagami-m fading severity
order = 50;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_int_blocks = zeros(length(gam ),length(U), length(m), length(famatype));
pout_calc_gl_blocks = zeros(length(gam_sim ),length(U), length(m), length(famatype));
pout_sim_blocks = zeros(length(gam_sim ),length(U), length(m), length(famatype));


%-----------------------------------------------------------------
%  Block correlation
%----------------------------------------------------------------
% Correlation matriz
Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));
% Eigenvalues
rho = sort(eig(Sigma_jakes),'descend');
% Correlation coefficent per block
deltab = 0.97;
% Number of domminant eigenvalues
Num_eig = sum(rho > N/100);

% Algorithm 1. L: vector with block sizes (Lb)
L = BlockCorrelation(N, rho, Num_eig, deltab);


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter =  length(U)  *  length(m) * length(famatype);
kiter = 0;

%-----------------------------------------------------------------
% Loop over  m
%-----------------------------------------------------------------
for km = 1:length(m)

    %-----------------------------------------------------------------
    % Loop over  U
    %-----------------------------------------------------------------
    for ku = 1:length(U)

        %-----------------------------------------------------------------
        % Loop over famatype 
        %-----------------------------------------------------------------
        for kfama = 1:length(famatype )

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);
     
            % Integral calculation
             pout_calc_int_blocks(:,ku,km,kfama) = CalcOutageFAMA(gam, L, deltab, U(ku), 'Integral', order, m(km), famatype(kfama), flag2int,flagmf, Inf);

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl_blocks(:,ku,km,kfama) = CalcOutageFAMA(gam_sim, L, deltab, U(ku), 'Quadrature', order, m(km),famatype(kfama),flag2int,flagmf, Inf);

            %-------------------------------------------------------------
            % Monte Carlo simulation
            % Loop over gam_sim
            %-------------------------------------------------------------

            for kg = 1:length (gam_sim)
                % Pout calculated
                pout = pout_calc_gl_blocks(kg,ku,km,kfama);

                % Avoid simulation for pout < 1e-5, as computation time would be very high
                if pout > 1e-5
                    NsamplesAdjusted = max( round( factor_sim * ( 1 / pout ) ), min_samples );
                    pout_sim_blocks(kg,ku,km,kfama) = SimOutage_BlocksFAMA(NsamplesAdjusted, gam_sim(kg), U(ku), sqrt(deltab), L, m(km), famatype(kfama),flagsimul, Inf);
                end
            end
        end
    end   
end
exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

Fig4_2curves
