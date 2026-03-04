%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots Multipexig gain vs Number of users for a network with:
% - cell-free with maximum ratio transmission (MRT) precoding 
% at base station and users with fluid antennas system (FAS)
% - fast and slow fluid antenna multiple access (FAMA)
% under Nakagami-m fading channels
% - spatial block-correlation model.
% 
% Computes: 
% - OP based on SIR using Gauss-Laguerre quadrature; 
% - Multiplexing gain for the OP calculated
% 
% Makes Figure 2.
%
% Parameters:
% N = 100; antBS = 2; m = 2; gamdB = 3; W = 1; d0 = 40; d = 100; alphaCF = 3;
% U = [4,5,6,7,8,9,10,15,20,25,30]
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

famatype = ["CFslow";"CFfast"];     %  FAMA type
N = 100;                   % Number of ports
W = 1;                      % Antenna size (wavelength normalized)

% Number of users
U = [4,5,6,7,8,9,10,15,20,25,30]; 


gamdB = 3;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)


m = 2;                    % Nakagami-m fading severity
order = 30;               % Order of GL quadrature

alphaCF = 3;
d0 = 40;  % obs.: 35 ficou bom
d = 100;
antBS = 2;

% flag, approx. mf = 1
flagmf = true;

% flag, 2nd integral
flag2int = true;


%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

mult_gain = zeros( size(U,1), size(U,2) ); 
pout = zeros( size(U,1), size(U,2) ); 


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

num_iter = length(U);
kiter = 0;


%-----------------------------------------------------------------
% Loop over number of Users
%-----------------------------------------------------------------
for ku = 1: length(U)

    % User feedback
    kiter = kiter + 1;
    disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);


    % OP - calculation - slow
    kfama = 1;
    orderSINR = order;
    pout(kfama, ku) = CalcOutageFAMA(gam, L, deltab, U (ku) , 'Quadrature', order, m, famatype(kfama),flag2int, flagmf, Inf, orderSINR, d0, d, alphaCF, antBS);

    % Multiplexing gain
    mult_gain(kfama, ku) = U(ku) * (1 - pout(kfama,ku));


    % OP - calculation - fast
    kfama = 2;
    orderSINR = 10 * order;
    pout(kfama, ku) = CalcOutageFAMA(gam, L, deltab, U (ku) , 'Quadrature', order, m, famatype(kfama),flag2int,flagmf, Inf, orderSINR, d0, d, alphaCF, antBS);

    % Multiplexing gain
    mult_gain(kfama, ku) = U(ku) * (1 - pout(kfama,ku));

end

exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

Figure2_curves
