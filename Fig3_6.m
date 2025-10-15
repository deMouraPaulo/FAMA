%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the Multiplexing Gain (Gm) vs Number of Users (U) for
% slow FAMA networks under Nakagami-m fading channels.
%
% Computes: OP based-SIR, integral expression.
% 
% Makes Figure 3.6
% Parameters: m = [1,2,3]; N = [200,300,400]; W = 2; gamma = 6 dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------

clc
clear
close all
addpath('Core/')


%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------

famatype = 'Slow';     % Set FAMA type

N = [200,300,400];     % Number of ports
W = 2;                 % Antenna size (wavelength normalized)
U = [2,3,4,5,6,7,8];   % Number of users


gamdB = 6;              % SINR threshold (dB)
gam = db2pow(gamdB);    % SINR threshold (linear scale)


m = [1,2,3];              % Nakagami-m fading severity
order = 30;               % Order of GL quadrature


%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------

mult_gain = zeros(length(N), length(m), length(U));
uoptimal = zeros (length(N), length(m));
mult_gain_opt = zeros (length(N), length(m));

%-----------------------------------------------------------------
% Multiplexing Gain calculation
%----------------------------------------------------------------- 

num_iter = length(N) * length(U) * length(m );
kiter = 0;
% Loop over  N
for kn = 1:length(N)

    % User feedback
    disp(['Iter "N": ' num2str(kn) ' out of ' num2str(length(N))]);

    % Average correlation coefficient
    dm = DeltaMed (N(kn), W);

    % Loop over  m
    for km = 1:length(m)

        % User feedback
        disp(['Iter "m": ' num2str(km) ' out of ' num2str(length(m))]);

        % Select equation to find the optimal number of users
        if m(km) == 1
            uopt = @(u,g,w,n)   u -  log(n)/( log( g/(1-1/(pi*w) )) )  - 1;
        elseif m(km) == 2
            uopt = @(u,g,w,n)   u - (log(2*u - 1) + log(n))/( 2*log( g/(1-1/(pi*w) )) )  - 1;
        elseif m(km) == 3
            uopt = @(u,g,w,n)   u - (log(4.5*(u^2- u) + 1) + log(n))/( 3*log( g/(1-1/(pi*w) )) )  - 1;
        end
        n = N(kn);
        g =  gam;
        w = W;
        fun = @(u) uopt(u,g,w,n);

        % Find root of the function "fun"
        uoptimal_dec  =  fzero(fun, 3);
        % Use the floor function
        uoptimal(kn,km) = floor(uoptimal_dec);
            

        % Loop over U
        for ku = 1:length(U)

            % User feedback
            kiter = kiter + 1;
            disp(['Total Iter : ' num2str(kiter) ' out of ' num2str(num_iter)]);

            % OP-based SIR
            pout = CalcOutageFAMA(gam, N(kn), dm, U(ku), 'Integral', order, m(km),famatype, 'SIR', Inf);
           
            % Multiplexing Gain
            mult_gain(kn,km,ku) = U(ku) * (1 - pout);

            % Multiplexing Gain for u_optimal
            if U(ku) == uoptimal(kn,km)
                mult_gain_opt(kn,km) = mult_gain(kn,km,ku);
            end
        end

    end

end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

% Colors for each "N" curve: "light blue", "light orange and "red"";
colors = [ [71, 147, 175]; [255, 196, 112]; [221, 87, 70] ]/255;

curveSIR = zeros(length(N));
for kn = 1:length(N)
    for km = 1:length(m)
        curveSIR(kn) = plot(U, squeeze(mult_gain(kn,km,:)),'s-' , ...
            'MarkerSize',10, 'Color',colors(kn,:), 'Linewidth', 2);
        hold on;
    end
end
curveSIR_opt = plot( uoptimal, mult_gain_opt, 'p', 'MarkerSize',16, ...
    'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6] );


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

% Colors for each "m" curve: "light blue", "light orange and "red"";
colors = [ [71, 147, 175]; [255, 196, 112]; [221, 87, 70] ]/255;

curveSIR = zeros(length(N));
for kn = 1:length(N)
    for km = 1:length(m)
        curveSIR(kn) = plot(U, squeeze(mult_gain(kn,km,:)),'s-' , ...
            'MarkerSize',10, 'Color',colors(kn,:), 'Linewidth', 2);
        hold on;
    end
end
curveSIR_opt = plot( uoptimal, mult_gain_opt, 'p', 'MarkerSize',18, ...
    'MarkerEdgeColor','r', 'MarkerFaceColor',[1 .6 .6],'LineWidth',2 );


%-------------------------------------------------------------------------
% Set figure parameters
%------------------------------------------------------------------------

grid on;
xlabel( 'Number of Users - $U$', 'Interpreter', 'Latex', 'FontSize', 24 );
ylabel( 'Multiplexing Gain - $\mathcal{G}_m$', 'Interpreter', 'Latex', 'FontSize', 24 );
xticks([2,4,6,8])
ylim( [0, 4] );
yticks([0,1,2,3,4])
set(gca, 'TickLabelInterpreter', 'latex','FontSize',24) 

legend([curveSIR(1), curveSIR(2), curveSIR(3),curveSIR_opt(1,1)],...
    {['$N = $ ',num2str(N(1))], ['$N = $ ',num2str(N(2))], ['$N = $ ',num2str(N(3))], '$\mathcal{G}_m(U^{\star})$'},...
    'Interpreter', 'Latex', 'FontSize', 24,'NumColumns',1, Location='northeast');

% Create textbox: parameters
annotation(figure(1),'textbox',...
    [0.78 0.52 0.10 0.10],...   
    'String',['$\gamma= $ ',num2str(gamdB),' dB', newline,'$W=$ ',num2str(W)],...
    'Interpreter','latex','FontSize',24,'HorizontalAlignment','center');

% Create textbox m = 1
annotation(figure(1),'textbox',...
    [0.446,0.783,0.0846,0.0373],...
    'String','$m = 1$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
annotation(figure(1),'ellipse',...
    [0.414,0.6006,0.021,0.219]);

% Create textbox m = 2
annotation(figure(1),'textbox',...
    [0.316,0.535,0.0846,0.03733],...
    'String','$m = 2$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
annotation(figure(1),'ellipse',...
    [0.2906,0.3951,0.02028,0.166]);

% Create textbox m = 3
annotation(figure(1),'textbox',...
    [0.1487,0.2216,0.0846,0.0373],...
    'String',{'$m = 3$'},...
    'Interpreter','latex',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
annotation(figure(1),'ellipse',...
    [0.2123,0.2675,0.0203,0.15]);


