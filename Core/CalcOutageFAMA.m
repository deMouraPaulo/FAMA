%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = CalcOutageFAMA(gamma_v, L, rho, U, method, order,... 
%     m, famatype, d0, d, alphaCF, Nant)
%
% This code serves as a supplementary tool  for theoretical analysis of
% "Fast, Slow and Opportunistic FAMA under Nakagami-m Fading Channels and
% Its Applications to Cell-Free Networks".
%
% Calculates theoretical OP for FAMA under Nakagami-m fading channels in
% different configurations:
% - systems: slow-FAMA, fast-FAMA, cell-free with MRT precoding for fast and
% slow FAMA;
% - correlation: constant correlation, spatial-block correlation;
% - categories of interference: interference-limited, noise-limited 
% (no interference), combined interference and noise;
% - methods: direct numerical integration, Gauss-Laguerre quadrature
% approximation (faster computation).
%
% ------------ Input Parameters -------------------------------------
%
% Specify FAMA type by setting "famatype":
%        famatype == 'Fast' -> fast-FAMA
%        famatype == 'Slow' -> slow-FAMA
%        famatype == 'CFfast' -> cell-free with MRT precoding fast-FAMA
%        famatype == 'CFslow' -> cell-free with MRT precoding slow-FAMA
%
% For constant correlation, set:
%           L = N, number of ports (scalar);
%           rho = delta, average power correlation coefficient.
%
% If the number of users is set to 1 (U == 1), then a single user 
% fluid antenna system is considered, and the OP based on SNR is computed. 
%
% If opbasedSINR == 'SINR' (and U ~= 1), OP based on SINR is computed, 
% In this case, the average received SNR (gamma_avg) must be given, 
% as noise is cosidered, additionaly to interference. 
%
% If "opbasedSINR" ~= 'SINR' (and U ~= 1), OP base on SIR is computed.
% 
% General parameters:
% - gamma_v: vector containing the SIR, SINR or SNR thresholds in linear scale; 
% - gamma_avg: scalar, average received SNR (used for OP based on SNR and SINR);
% - L: vector containing the block sizes of the correlation approximation;      
% - rho: power correlation coefficient used in the block-diagonal, 
%        approximation (suitable range: 0.95 to 0.99);
% - U: number of users (scalar);
% - method: string that specifies the theoretical expression used to compute the OP 
%      If method == 'Integral', then the double or triple integral expressions 
%           are used for direct numerical integration;
%      If method == 'Quadrature', the integrals are solved by Gauss-Laguerre
%           quadrature;
% - order: order of the quadrature approximation 
%          (only if method == 'Quadrature' );
% - orderSINR : order of the 1st quadrature approximation 
%           (only if method == 'Quadrature' AND opbasedSINR == 'SINR');
% - m: Nakagami-m fading severity (integer);
%
%  Parameters for free-cell FAMA with MRT precoding:
%       d0 - distance from serving BS to the target user;
%       d - distance from interfering BSs to the target use, 
%               (d is the same for all interfering BS);
%       alphaCF - path loss exponent;
%       Nant - number of antennas in a given BS.
%
%------------------- Output -------------------------------------
%
% - pout: vector containing the OP (same size as gamma_v), i.e., 
%         pout(k) = P(SIR < gamma_v(k)).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = CalcOutageFAMA(gamma_v, L, rho, U, method, order, m, famatype,...
    opbasedSINR, gamma_avg, orderSINR, d0, d, alphaCF, Nant)

% -------------- A type of FAMA must be indicated -----------------
% ---------------------------------------------------------------------
if strcmp(famatype, 'Fast') || strcmp(famatype, 'Slow') || strcmp(famatype, 'CFfast') || strcmp(famatype, 'CFslow')

    %--------------- Single user case ------------------------------------
    %---------------------------------------------------------------------
    if (U == 1 )  % FAS OP 
        % Pre-allocating
        pout = ones(1,length(gamma_v));

        if strcmp(famatype, 'CFfast') || strcmp(famatype, 'CFslow')
             m = m * Nant;
        end
        gamma_ratio = gamma_v ./ gamma_avg;

        % loop over blocks
        for kb = 1:length(L) 

            fun = @(rb) rb.^(m-1) .* exp(-rb/2).*...
                (1-marcumq(sqrt(rho*rb./(1-rho)),sqrt( (2*m*gamma_ratio)./(1-rho) ), m)).^L(kb);

            pout = pout.* (1/(2^m * gamma(m)) * integral(fun, 0, inf, 'ArrayValued',true));

        end
        

    %--------------- Multi user case -------------------------------------
    %---------------------------------------------------------------------
    else
        if strcmp(famatype, 'Slow')
            Util = m * (U-1);   

        elseif strcmp (famatype, 'CFslow')
            Util = m * (U-1);
            gamma_v = gamma_v .* (( d0 / d ).^alphaCF);
            m = m * Nant;

        elseif strcmp(famatype, 'Fast')
            Util = 1;  % approximation for OP-SIR
            gamma_v = gamma_v .* ( m * (U - 2) + 1);

        elseif  strcmp(famatype, 'CFfast')
            Util = 1;  % approximation for OP-SIR
            gamma_v =  gamma_v .* ((m*(U-2) + 1) * ( d0 / d ).^alphaCF );
            m = m * Nant;
        end

        % ----------------- SIR based OP - 2 integrals -----------------
        % --------------------------------------------------------------
        if ~strcmp(opbasedSINR, 'SINR')

            %---------------- Direct integration method  -------------
            if strcmp(method, 'Integral')
                % Pre-allocating
                pout = ones(1,length(gamma_v));

                % Loop over SIR thresholds (gamma_v)
                for kg = 1:length(pout)

                    % Loop over blocks
                    for kmu = 1:length(L)

                        % Define integrand with particularized values
                        fint = @(r_b, rtil_b) r_b.^(m-1).* rtil_b.^(Util-1).*exp(-(r_b+rtil_b)/2)...
                            .*Gfun(r_b, rtil_b, gamma_v(kg), Util, rho, m).^L(kmu);
                   
                        % Perform numerical double integration
                        pout(kg) = pout(kg) * ...
                            (integral2(fint, 0, inf, 0, inf))/(2^(m+Util)*gamma(m)*gamma(Util));

                    end

                end  % END loop over SIR thresholds (gamma_v) 
                % END method, 'Integral'

            %---------------- Quadrature approximation  -----------------
            elseif strcmp(method, 'Quadrature')

                % The function "Gauss-Laguerre" is used to avoid complex roots
                [x, w] = GaussLaguerre(order, m-1);
                [xtil, wtil] = GaussLaguerre(order, Util-1);

                % Store each summand in a 3D tensor
                Gmatrix = zeros(order, order, length(gamma_v));
                for ki = 1:order
                    for kj = 1:order
                        Gmatrix(ki,kj,:) = Gfun(2*x(kj), 2*xtil(ki), gamma_v, Util, rho, m);
                    end
                end

                % Computes outage probability
                pout = ones(length(gamma_v),1);
                for kmu = 1:length(L)
                    pout = pout.* squeeze(sum(wtil.*(Gmatrix.^L(kmu)).*w.',[1 2]))/(gamma(m)*gamma(Util));
                end

            else
                error('Integration method must be indicated')
            end % integration method
        
        %---------------- OP based SINR: 3 integrals ------------------
        %-------------------------------------------------------------
        elseif strcmp(opbasedSINR, 'SINR')

            % calculates: SINR / SINR average
            if strcmp(famatype, 'Slow')
                gamma_ratio = gamma_v ./ gamma_avg;
            elseif strcmp (famatype, 'CFslow')
                gamma_ratio = gamma_v ./ (gamma_avg .* (( d0 / d ).^alphaCF));
            elseif strcmp(famatype, 'Fast')
                gamma_ratio = gamma_v ./ (gamma_avg * ( m * (U - 2) + 1));
            elseif strcmp(famatype, 'CFfast')
                gamma_ratio = gamma_v ./ (gamma_avg * ( m * (U - 2) + 1) * ( d0 / d ).^alphaCF );
            end

            % Util = mtil for f-FAMA
            if strcmp(famatype, 'Fast') ||  strcmp(famatype, 'CFfast')
                Util = (m * (U-1)) / (m * (U-2) + 1);
            end


            %----------- Direct integration method  -------------------
            if strcmp(method, 'Integral')

                % 1st integration uses Matlab function "integral"
                %
                % 2nd and 3th integrations use Matlab function "trapz", with
                % these parameters for numerical integration: 
                points_integral = 100;
                rb_max = 15*m;
                r_b = linspace( 0, rb_max, points_integral );
                rtil_b_max = 15*Util;
                % lower limit is not zero to avoid division by zero
                rtil_b = linspace( 1e-8, rtil_b_max, points_integral );

               % constant
                df = rho / ( 1 - rho );
                
                % Pre-allocating
                pout = ones(1,length(gamma_v));

                % Loop over SINR thresholds (gamma_v)
                for kg = 1:length(pout)

                    % Loop over blocks
                    for kmu = 1:length(L)       

                        % 3th integrand
                        fr_b = zeros( length( r_b ), 1 );
                        for ib = 1 : length( r_b )

                            % 2nd integrand
                            frtil_b = zeros( length( rtil_b ), 1 );
                            for jb = 1 : length( rtil_b )

                                % 1st function to integrate
                                fint = @(yk) ( yk / ( ( df * rtil_b( jb ) ) ) ).^( 0.5 * ( Util - 1 ) ) .*....
                                    exp( sqrt( df  * rtil_b( jb ) * yk ) - 0.5 *( yk +  df * rtil_b( jb )  ) ) .*...
                                    besseli( Util - 1, sqrt( df  * rtil_b( jb ) * yk ), 1) .*...
                                    marcumq( sqrt( df * r_b( ib ) ), ...
                                    sqrt( gamma_v(kg) * yk  +  ( 2 * m * gamma_ratio(kg) ) / ( 1 - rho )  ), m );

                                % Perform 1st numerical integration
                                f1 = ( 1 - 0.5 * integral(fint, 0, Inf) ).^( L(kmu) );

                                % 1st integral * exponential and power terms
                                frtil_b( jb, 1 ) = f1 * exp( - 0.5 * rtil_b( jb ) ) * rtil_b( jb )^( Util - 1 );
                            end
                            % calculates 2nd integral
                            f2 = trapz( rtil_b, frtil_b );

                            % 2nd integral * exponential and power terms
                            fr_b( ib, 1 ) = f2 * exp( - 0.5 * r_b( ib ) ) * r_b( ib )^( m - 1 );
                        end

                        % 3nd integral -> OP
                        pout(kg) = pout(kg) .*( trapz( r_b, fr_b ) / (2^(Util+m) * gamma(m) * gamma(Util)) );

                    end % end Loop over blocks (Lb)
                end % end loop over SINR thresholds (gamma_v)

            %--------- Quadrature approximation for OP-SINR -------------
            %------------------------------------------------------------
            elseif strcmp(method, 'Quadrature') 
                %
                % set order of GLQ approximation
                n_i = orderSINR; % for f-FAMA n_i should be higher than "order" due to MarcumQ term
                n_j = order;   
                n_l = order;

                % The function "Gauss-Laguerre" is used to avoid complex roots
            	[a_i, w_i] = GaussLaguerre( n_i, 0 );
               	[a_j, w_j] = GaussLaguerre( n_j, Util - 1 );
                [a_l, w_l] = GaussLaguerre( n_l, m - 1 );                
                
                % constant
                df = rho / ( 1 - rho );

                % Pre-allocating
                pout = ones(1,length(gamma_v));

                % Loop over SINR thresholds (gamma_v)
                for kg = 1:length(pout)

                    % Loop over blocks
                    for kmu = 1:length(L)

                        f3 = zeros(n_l, 1);
                        for l = 1 : n_l   % loop over 3th GLQ (n_l)

                            % Marcum-Q Function
                            marc_t = marcumq( sqrt( 2 * df * a_l( l ) ),...
                                sqrt( (2*gamma_v(kg) * a_i) + ( 2 * m * gamma_ratio(kg) ) / ( 1 - rho )  ), m );

                            f2 = zeros( n_j, 1 );
                            for jj = 1 : n_j  % loop over 2nd GLQ (n_j)

                                % 1st integrand
                                fint = ( a_i / ( df.* a_j( jj ) ) ).^( 0.5 * ( Util - 1 ) ).*...
                                    exp(sqrt( 4 * df * a_i * a_j( jj ) ) -  df.* a_j( jj )  ) .*...
                                    besseli( Util - 1, sqrt( 4 * df * a_i * a_j( jj ) ), 1 ).*...
                                    marc_t;

                                % 1st GLQ
                                f1 = ( 1 - sum( fint .* w_i ) ).^( L(kmu) );

                                % 1st GLQ * 2nd weight
                                f2( jj, 1 ) = f1 * w_j( jj );
                          
                            end % End loop over 2nd GLQ (n_j)

                            % 2nd GLQ * 3th weight
                            f3( l, 1 ) = sum( f2 )  * w_l( l );

                        end % End loop over 3th GLQ (n_l)

                        % 3th GLQ -> OP
                        pout(kg) = pout(kg) .*  (sum( f3 ) / ( gamma( m ) * gamma( Util ) )) ;

                    end  % end Loop over blocks (Lb)
                end % end loop over SINR thresholds (gamma_v)

            else
                error('Integration method must be indicated')   
            end % if method - integral or quadrature
            
        end  % END opbased 'SINR' or 'SIR'
    end % END single (U == 1) or multi-user 


else
    error('Type of FAMA must be indicaded: Fast, Slow, CFfast, CFslow')
end

    %---------------------------------------------------------------------
    % Auxiliary function to compute G 
    %---------------------------------------------------------------------
    function [G] = Gfun(r_b, rtil_b, gam, Util, rho, m)

        summatory = 0;

        % To avoid overflow or loss of accuray, the modified Bessel
        % function of the first kind, besseli(nu,Z,scale), was scaled by
        % the factor exp(-abs(real(Z))), with scale = 1. 
        % Therefore, the exponential term was multiplied by exp(+Z).
        %
        % gammaln(A) = log(gamma(A)) 
        % The gammaln command avoids the underflow and overflow that 
        % may occur if it is computed directly using log(gamma(A)).

        for k = 0:(m+Util-2)
            for j = 0:(m+Util-k-2)
                summatory = summatory + exp(gammaln(m+Util-k-1) - gammaln(m+Util-j-k-1)...
                    - gammaln(j+1) + ((j+k)/2)*log(r_b./rtil_b)) .* (gam+1).^k ...
                    .* gam.^((j-k)/2) .* ...
                    besseli(1-m+j+k, rho*sqrt(gam.*r_b.*rtil_b)/(1-rho)./(gam+1),1)...
                  .* exp((rho/(1-rho)./(gam+1)).*...
                        (sqrt(gam.*r_b.*rtil_b)-(gam.*rtil_b+r_b)/2));
            end
        end

        G = marcumq(sqrt(rho*gam.*rtil_b/(1-rho)./(gam+1)), ...
                       sqrt(rho*r_b/(1-rho)./(gam+1)),Util) - ...
               (1./(gam+1)).^(m+Util-1).*((r_b./(gam.*rtil_b)).^((1-m)/2)).*summatory;

        G(isnan(G)) = 0; 

    end

    %---------------------------------------------------------------------
    % Gauss-Laguerre quadrature
    %---------------------------------------------------------------------
    function [x, w] = GaussLaguerre(n, alpha)
        % This function determines the abscisas (x) and weights (w) for the
        % Gauss-Laguerre quadrature of order n>1, on the interval [0, +infinity].
            % Unlike the function 'GaussLaguerre', this function is valid for
            % n>=34. This is due to the fact that the companion matrix (of the n'th
            % degree Laguerre polynomial) is now constructed as a symmetrical
            % matrix, guaranteeing that all the eigenvalues (roots) will be real.
                       
        % © Geert Van Damme
        % geert@vandamme-iliano.be
        % February 21, 2010    
        % Building the companion matrix CM
            % CM is such that det(xI-CM)=L_n(x), with L_n the Laguerree polynomial
            % under consideration. Moreover, CM will be constructed in such a way
            % that it is symmetrical.
        i   = 1:n;
        a   = (2*i-1) + alpha;
        b   = sqrt( i(1:n-1) .* ((1:n-1) + alpha) );
        CM  = diag(a) + diag(b,1) + diag(b,-1);
        % Determining the abscissas (x) and weights (w)
            % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
            %   characteristic polynomial, i.d. the eigenvalues of CM;
            % - the weights can be derived from the corresponding eigenvectors.
        [V, Lg]   = eig(CM);
        [x, ind] = sort(diag(Lg));
        V       = V(:,ind)';
        w       = gamma(alpha+1) .* V(:,1).^2;
        
    end

end


