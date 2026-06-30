%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the Multipexig gain vs Number of users.
% Fast, Slow and Opportunistic FAMA networks with block-correlation model.
% 
% Makes Figure 4.3
%
%
% Parameters:
% W = 1, N = 100, gamma = 3dB, m = 2
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

% Number of users - slow FAMA ; Number of users - fast FAMA
U = round ([linspace(2,15,14); linspace(10,150,14)]); % (famatype, U) 

% Pool of users (famatype, U,  M)
Mpool = zeros(size(U,1),size(U,2),3);
Mpool(:,:,1) = round (1.5 * U);
Mpool(:,:,2) = round (2.0 * U);
Mpool(:,:,3) = round (2.5 * U);

% Transpose U (famatype, U )
% U = U';


    
gamdB = -3;              % SIR threshold (dB)
gam = db2pow(gamdB);    % SIR threshold (linear scale)


m = 2;                    % Nakagami-m fading severity
order = 30;               % Order of GL quadrature



%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

mult_gain = zeros( size(U,1), size(U,2) ); 
mult_gain_OFAMA = zeros( size(U,1), size(U,2), size(Mpool,3) );


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

num_iter = numel(Mpool);
kiter = 0;


%-----------------------------------------------------------------
% Loop over  famatype
%-----------------------------------------------------------------
for kfama = 1:length(famatype)

    %-----------------------------------------------------------------
    % Loop over number of Users
    %-----------------------------------------------------------------
    for ku = 1: size (U,2) 
     
        % OP - calculation
        pout = CalcOutageFAMA(gam, L, deltab, U (kfama, ku) , 'Quadrature', order, m, famatype(kfama),flag2int,flagmf, Inf);

        % Multiplexing gain: s-FAMA and f-FAMA
        mult_gain(kfama, ku) = U(kfama, ku) * (1 - pout);

        %-----------------------------------------------------------------
        % Multiplexing gain O-FAMA
        % Loop over  Mpool
        %-----------------------------------------------------------------
        for kpool = 1:size(Mpool,3)

            % User feedback
            kiter = kiter + 1;
            disp(['Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            mult_gain_OFAMA(kfama, ku, kpool) = MultiplexGain(Mpool(kfama, ku, kpool), U(kfama, ku), pout );

        end

    end
end
exec_time = toc/60;
disp (['Execution time: ', num2str(exec_time), ' min'])


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 

Fig4_3curves
