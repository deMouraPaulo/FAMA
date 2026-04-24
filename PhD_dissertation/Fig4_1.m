%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs Number of ports.
% Fast and slow FAMA networks, constant and block-correlation models.
% 
% Makes Figure 4.1
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
% Parameters
%-------------------------------------------------------------------------

famatype = ["Slow";"Fast"];     %  FAMA type
N = 5:5:200;                   % Number of ports
Nsimul = [5,10,20,60,100,140,180,200];    % Number of ports for simulation
W = [1,3];                      % Antenna size (wavelength normalized)
U = 5;                          % Number of users

min_samples = 5e5;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e2;            % Multiplication factor for simulation 
     
gamdB = -3;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)


m = 2;                    % Nakagami-m fading severity
order = 30;               % Order of GL quadrature
order3rd = order;         % order of 3rd sum for GLQ calculation


flagcalc = true;  % (ffama calc, 2 integrals)
flagsim  = true;  % (ffama sim, approximated, Gaussian RVs)

flagmf = true;  % (m_fast = 1)



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_calc_int = zeros(length(N), length(W), length(famatype));
pout_calc_int_blocks = zeros(length(N), length(W), length(famatype));

pout_calc_gl = zeros(length(Nsimul), length(W), length(famatype));
pout_sim = zeros(length(Nsimul), length(W), length(famatype));

pout_calc_gl_blocks = zeros(length(Nsimul), length(W), length(famatype));
pout_sim_blocks = zeros(length(Nsimul), length(W), length(famatype));

pout_sim_jakes = zeros(length(Nsimul), length(W), length(famatype));


%-----------------------------------------------------------------
% Outage Probabilities calculation
%----------------------------------------------------------------- 

num_iter = length(N) *  length(W) * length(famatype);
kiter = 0;
% Loop over  W
for kw = 1:length(W)

    % User feedback
    disp(['Iter "W": ' num2str(kw) ' out of ' num2str(length(W))]);

    % Loop over  famatype
    for kfama = 1:length(famatype)

        %-----------------------------------------------------------------
        % Loop over N - Integral calculation
        %-----------------------------------------------------------------
        for kn = 1:length(N )

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            %-----------------------------------------------------------------
            % Constant Correlation
            %-----------------------------------------------------------------
            % Average correlation coefficient
            dm = DeltaMed (N(kn), W(kw));

            % Integral calculation
            pout_calc_int(kn,kw,kfama) = CalcOutageFAMA(gam, N(kn), dm, U, 'Integral', order, m, famatype(kfama), flagcalc, flagmf, Inf);        
            
            %-----------------------------------------------------------------
            %  Block correlation
            %----------------------------------------------------------------
            % Correlation matriz
            Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N(kn)-1)*W(kw)/(N(kn)-1)));
            % Eigenvalues 
            rho = sort(eig(Sigma_jakes),'descend');
            % Correlation coefficent per block
            deltab = 0.97;
            % Number of domminant eigenvalues
            Num_eig = sum(rho > N(kn)/100);
            % Algorithm 1. L: vector with block sizes (Lb)
            L = BlockCorrelation(N(kn), rho, Num_eig, deltab);

            % Integral calculation
            pout_calc_int_blocks(kn,kw,kfama) = CalcOutageFAMA(gam, L, deltab, U, 'Quadrature', order, m, famatype(kfama), flagcalc, flagmf, Inf,order3rd);

        end

        %-----------------------------------------------------------------
        % Loop over Nsimul - Simulation and Quadrature Calculation
        %-----------------------------------------------------------------
        for kn = 1:length( Nsimul )

            %-----------------------------------------------------------------
            % Constant Correlation
            %-----------------------------------------------------------------
            % Average correlation coefficient
            dm = DeltaMed (Nsimul(kn), W(kw));

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl(kn,kw,kfama) = CalcOutageFAMA(gam, Nsimul(kn), dm, U, 'Quadrature', order, m,famatype(kfama), flagcalc, flagmf, Inf,order3rd);
            
            %------------ Monte Carlo simulation -------------------------

            % Pout calculated
             pout = pout_calc_gl(kn,kw,kfama);

            % Avoid simulation for pout < 0.9e-5, as computation time would be very high
            if pout > 0.9e-5
                NsamplesAdjusted = max( round( factor_sim * ( 1 / pout ) ), min_samples );
                pout_sim(kn,kw,kfama) = SimOutage_BlocksFAMA(NsamplesAdjusted, gam, U, sqrt(dm), Nsimul(kn), m,famatype(kfama), flagsim, Inf);
            end
            %-----------------------------------------------------------------
            %  Block correlation
            %----------------------------------------------------------------
            % Correlation matriz
            Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:Nsimul(kn)-1)*W(kw)/(Nsimul(kn)-1)));
            % Eigenvalues 
            rho = sort(eig(Sigma_jakes),'descend');
            % Correlation coefficent per block
            deltab = 0.97;
            % Number of domminant eigenvalues
            Num_eig = sum(rho > Nsimul(kn)/100);
            % Algorithm 1. L: vector with block sizes (Lb)
            L = BlockCorrelation(Nsimul(kn), rho, Num_eig, deltab);

            % Gauss-Laguerre Quadrature calculation
            pout_calc_gl_blocks(kn,kw,kfama) = CalcOutageFAMA(gam, L, deltab, U, 'Quadrature', order, m,famatype(kfama), flagcalc, flagmf, Inf,order3rd);

            % Monte Carlo simulation
            NsamplesAdjusted = max( round( factor_sim * ( 1 / pout_calc_gl_blocks(kn,kw,kfama) ) ), min_samples );
            % based on SIR
            pout_sim_blocks(kn,kw,kfama) = SimOutage_BlocksFAMA(NsamplesAdjusted, gam, U, sqrt(deltab), L, m, famatype(kfama),flagsim, Inf);
            % based on Jake's matriz
            pout_sim_jakes(kn,kw,kfama) = SimOutageFAMA(NsamplesAdjusted, gam, Sigma_jakes, U, m, famatype(kfama),Inf);
        end
    end   
end
exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
Fig4_1curves


