%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Typical fading for FAS
% Plots channel magnitude x linear distance
%
% Makes Figure 2.1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization
close all
clear

% Parameters
users = 1; % number of users
N = 100; % number of ports
W = 2;  % antenna size

% Correlation matrix, Jake's model
Sigma = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));

% Eigenvectors and eigenvalues
[V,Lambda] = eig(Sigma);

% Sort eigenvalues in descending order
[lambda, index] = sort(diag(Lambda),'descend');

% Only the eigenvalues larger than a small tolerance are stored
lambda = lambda(lambda > 1e-5); % 1e-5
 
% Re-arranges the corresponding eigenvectors
V = V(:,index(1:length(lambda)));

% Square-root of correlation matrix
A = V*diag(sqrt(lambda));

% Number of independent Gaussian RVs (rank of Sigma)
Ngen = length(lambda);

% Channel magnitude
% "like=1i" generrate complex Gaussian RV CN(0,1)
ch = abs(A * randn(Ngen, users, like=1i) )./sqrt(2); 

% Convert to dB
chdB = mag2db (ch);

% Distance from port 1, in waveleghts
k = 1:1:N;

% Plot channel samples
dist = (k-1)*W / (N-1); % dist from port 1, in wavelenghs
plot(dist, chdB, LineWidth=3);
hold on;

% Figure labels
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20)
ylabel("channel magnitude (dB)", 'FontSize', 20,'Interpreter','latex');
xlabel("linear distance (in wavelenghts)",'FontSize', 20, 'Interpreter','latex')