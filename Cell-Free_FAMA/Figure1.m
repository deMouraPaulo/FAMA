%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP vs SIR threshold for a network with:
% - cell-free with maximum ratio transmission (MRT) precoding 
% at base station and users with fluid antennas system (FAS);
% - fast and slow fluid antenna multiple access (FAMA) under Nakagami-m 
% fading channels;
% - constant correlation model.
%
% Computes: 
% - OP based on SIR using Gauss-Laguerre quadrature;
% - OP based on SIR using integral;
% - OP based on SIR, Monte Carlo simulation;
% 
% Makes Figure 1
%
% Parameters:
% N = 10; antBS = [2,8] m = 2; U = 4; W = 1; d0 = 20; d = 50; alphaCF = 3;
% SIR threshold: from 5 to 20 dB.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------
tic
clc
close all
clear

addpath(fullfile('..','Core'))


%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------


N = 10;                     % Number of ports

W = 1;                       % Antenna size (wavelength normalized)
U = 3+1;                   % Number of users

         
gamdB = linspace(5,20,10);        % SIR threshold (dB)
gam = db2pow(gamdB);               % SIR threshold (linear scale)


m = 2;                      % Nakagami-m fading serverity
order = 30;                    % Orderm of GL quadrature
orderSINR = 30;             % Order of 1st GL quadrature for OP-SINR


alphaCF = 3;
d0 = 20;
d = 50;
antBS = [2,8];

% flag, approx. mf = 1
flagmf = true;

% flag, 2nd integral
flag2int = true;

% flag, simulation, 
flagsimul = true; % (mf = 1, approx. method, Gaussian RVs)


min_samples = 5e5;          % Minimum number of samples for Monte-Carlo simulation
factor_sim = 1e2;            % Multiplication factor for simulation 

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

pout_simCFslow = zeros(length(antBS),length(gam));
pout_simCFfast = zeros(length(antBS),length(gam));
pout_calcfast_int = zeros(length(antBS),length(gam));
pout_calcslow_int = zeros(length(antBS),length(gam));
pout_calcfast_quad = zeros(length(antBS),length(gam));
pout_calcslow_quad = zeros(length(antBS),length(gam));

% Average correlation coefficient
dm = DeltaMed (N, W);


%-----------------------------------------------------------------
% Outage Probabilities calculation/simulation
%-----------------------------------------------------------------

num_iter =  length(antBS)  *  length(gam);
kiter = 0;

for ka = 1:length(antBS)

    famatype = "CFslow";
    % Calculation
    pout_calcslow_int(ka,:) = CalcOutageFAMA(gam, N, dm, U, 'Integral',order, m,famatype,flag2int,flagmf,  Inf , orderSINR, d0, d, alphaCF, antBS(ka));
    pout_calcslow_quad(ka,:) = CalcOutageFAMA(gam, N, dm, U, 'Quadrature',order, m,famatype,flag2int,flagmf,  Inf , orderSINR, d0, d, alphaCF, antBS(ka));


    % Simulation only for pout < 1e-6, to avoid high computation time
    for kg = 1:length(gam)

        % User feedback
        kiter = kiter + 1;
        disp(['Iter ' num2str(kiter) ' out of ' num2str(num_iter)]);

        pout = pout_calcslow_int(ka,kg);
        if pout > 1e-6
            Nsamples = max( round( factor_sim * ( 1 / pout ) ), min_samples );
            pout_simCFslow(ka,kg) = SimOutage_BlocksFAMA(Nsamples,gam(kg), U, sqrt(dm), N, m,famatype,flagsimul, Inf, d0, d, alphaCF,antBS(ka));
        end
    end


    famatype = "CFfast";
    % Calculation
    pout_calcfast_int(ka,:) = CalcOutageFAMA(gam, N, dm, U, 'Integral',order, m,famatype,flag2int,flagmf,Inf , orderSINR, d0, d, alphaCF, antBS(ka));
    pout_calcfast_quad(ka,:) = CalcOutageFAMA(gam, N, dm, U, 'Quadrature',order, m,famatype,flag2int,flagmf,  Inf , orderSINR, d0, d, alphaCF, antBS(ka));


    % Simulation only for pout < 1e-6, to avoid high computation time
    for kg = 1:length(gam)
        pout = pout_calcfast_int(ka,kg);
        if pout > 1e-6
            Nsamples = max( round( factor_sim * ( 1 / pout ) ), min_samples );
            pout_simCFfast(ka,kg) = SimOutage_BlocksFAMA(Nsamples,gam(kg), U, sqrt(dm), N, m,famatype, flagsimul, Inf, d0, d, alphaCF,antBS(ka));
        end
    end

end
exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

Figure1_curves

