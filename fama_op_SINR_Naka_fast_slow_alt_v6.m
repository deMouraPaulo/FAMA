function out_prob = fama_op_SINR_Naka_fast_slow_alt_v6( num_ports, num_users, gamma_th, delta, gamma_avg, m_nkg, points_integral,fast, plot_graf, ~ )
% FUNCTION NAME:
%   fama_op_SINR_Naka_fast_slow_alt_v6
%
% DESCRIPTION:
%   Calculates the outage probability based on the exact expression SINR
%   for slow and fast FAMA over Nakagami-m channels [new articles]
%   All channels has the same Nakagami-m fading factor m_nkg.
%
%   V6: Integral em yk simbólica para evitar ymax diferente,s dependendo dos
%   parâmetros, mas o tempo de integração é elevado.
%
%   v5: CORRECT calculation ERRORS: f_total ( isnan (f_total)) = 0
%   - Alterado ymax_inicial do f-FAMA de 2 para 5e-2, pois melhorar a
%   precisão, em especial, para gamma_th = 20 dB. Evita soma de valores
%   nulos
%   Alterado r_final, de 10*m_nkg para 15*m_nkg, m > 1, então ficou
%   r_final = 15*m_nkg para todos
%
%   V4: usa expressão alternativa para OP, com função de Bessel substituída
%   por uma integral. Isso ajuda a evitar problemas numéricos, pois a
%   função de Bessel tende infinito quando seu argumento tende ao infinito.
%   Assim um número muito grande multiplicado por outro pequena (na
%   integral) provoca problema numérico.
%   - Retirado parte do código (versão anterior) que avalia se a 
%   função de Bessel assume valores infintos, defido a nova expressão.
%   - Valores iniciais de rp, r, alterados de (1e-3, 1e-3) p/ (1e-8, 0)
%   - Criada função para determinar rp_max
%   - Criadas funções para plotagem de gráficos para analisar correação dos
%   cálculos
%   - Melhorados organização e comentários.
%   - Obs.: usar valores realistas de m_nkg e num_users, compatíveis com fast
%   e slow FAMA para os limites de integração ficarem adequados. Observar
%   gráficos para verificar os limites.
%   - check calculation errors
%
% INPUT:
%   num_ports - (int) Number of ports
%   num_users - (int) Number of users
%   gamma_th - (double) SINR threshold
%   delta - (double) Correlation factor
%   gamma_avg - Average SNR = ( 2 * m_nkg * sigma_s^2 * sigma^2 ) / ( sigma_n^2 );
%       sigma_g - (double) Square root of channel power
%       sigma_s - (double) Square root of symbol power
%       sigma_n - (double) Square root of noise power
%   m_nkg - (int) Nakagami-m fading intensity factor
%   points_integral - (int) Number of points used in the integrations
%   fast - (boolean), if fast-FAMA = TRUE, if slow-FAMA = FALSE 
%   plot_graf - (boolean), if plot_graf = TRUE plot grafs of integrals
%   new_expression - (boolean), NOT USED
%
% OUTPUT:
%   out_prob - (double) Outage probability
%


% Test: "fast FAMA" or "slow FAMA"
if fast  == true % fast FAMA
    % for fast-FAMA u_eq = m til, fading factor of overall interference
    u_eq =  1 / ( 1 - (m_nkg-1) / (m_nkg * (num_users -1)) );
    % for slow-FAMA, gamma_eq = gamma_th * U hat
    gamma_eq = gamma_th * ( m_nkg * (num_users - 2) + 1);
else  % "slow FAMA"
    % for slow-fama, u_eq = U til = sum of m_nkg for all users
    u_eq = (num_users - 1) * m_nkg;
    % for slow FAMA, gamma_eq = gamma_th
    gamma_eq = gamma_th;
end

% Constants
k1 = 1 / ( 2^( u_eq + m_nkg  ) * gamma( m_nkg ) * gamma( u_eq ) );
d_f = delta / ( 1 - delta );


% Range for varible "r". Desired user (last integral)
r_max = 15*m_nkg;
r = linspace( 0, r_max, points_integral );


% Range for varible "rp". Interference user
rpmax_inicial = 100;
thershold = 1e-5;
rp_max = find_maxrp ( u_eq, thershold, rpmax_inicial, points_integral );
% Final domain for "rp"
rp = linspace( 1e-8, rp_max, points_integral );
 % obs: initial value cannot be = 0, to avoid division by zero


% Symbolic integration of yk, so ymax can be large
ymax = 1000;


% 3th integrand
fr_aux = zeros( length( r ), 1 );
for i = 1 : length( r ) 
        
    % 2nd integrand
    frp_aux = zeros( length( rp ), 1 );
    for j = 1 : length( rp )  

       % function to integrate
        f_total = @(yk) ( yk / ( ( d_f.* rp( j ) ) ) ).^( 0.5 * ( u_eq - 1 ) ) .*....
           exp( -0.5 * ( yk + ( ( delta.* rp( j ) )./( 1 - delta ) ) ) ) .*...
           besseli( u_eq - 1, sqrt( yk * ( delta.* rp( j ) ) / ( 1 - delta ) ) ) .*...
           marcumq( sqrt( ( d_f * r( i ) ) ), ...
           sqrt( ( gamma_eq * yk ) + ( ( ( 2 * m_nkg * gamma_th ) / ( gamma_avg * ( 1 - delta ) ) ) ) ), m_nkg );

       % 1st integral
       intyk = integral(f_total, 0, ymax);

       % calculates (1 - 0.5 * 1st integral).^( num_ports )
       int_m_t = ( 1 - 0.5 * intyk ).^( num_ports );

       % 1st integral X exponential terms
       frp_aux( j, 1 ) = int_m_t * exp( - 0.5 * rp( j ) ) .* rp( j )^( u_eq - 1 );

    end % for j = 1 : length( rp ) - 2nd integrand
        
        
    % plot 2ª INT
    if plot_graf
        plot_2nd_int (i, r(i), length(r), rp, frp_aux)
    end

    % calculates 2nd integral
    Frp = trapz( rp, frp_aux );

    % 2nd integral X exponential terms
    fr_aux( i, 1 ) = Frp * exp( -r( i ) / 2 ) * r( i )^( m_nkg - 1 );

end % END for i = 1 : length( r ) - 3th integrand

    
% plot 3th integral 
if plot_graf
    plot_3th_int (r, fr_aux)
end

% OP - FINAL RESULT

% 3th integral
out_prob = k1 * trapz( r, fr_aux );

end % END MAIN function


function rp_max = find_maxrp ( u_eq, gl_th, rpmax_inicial, points_integral )
% DESCRIPTION:
%    Find the maximum value of "r" used in the integration.
%
% Obs: "f" is a suitable function evaluated to find r_max

    rp = linspace( 1e-8, rpmax_inicial, points_integral );
    f = exp( -rp / 2 ) .* rp .^( u_eq - 1 );
    
    % Find rp_max

    % f_total max
    [ ~, i_max] = max (f);
    % threshold for f_total
    f_th = gl_th * max( f );
    % Find min(index) from colums f_total(i_max) to f_total(end)
    [ ~, i_th ] = min( abs( f(i_max:end) - f_th ) );
    % r_max at the right of f_total(i_max) 
    % This is suitable for functions that increase and then decrease
    if i_th + i_max < length (rp)
        rp_max = rp ( i_th + i_max );
    else
        rp_max = rp ( length(r) );
    end
end % END function find_maxrp

function plot_2nd_int (i, ri,  lengthr, rp, frp_aux)
% Descrpition: plot 2nd integral
% vectors: rp, frp_aux;


if or( i==10, i==(lengthr-50) )
    plot (rp,frp_aux);
    str2 = {['2nd Integral (pos, valor): r = (' ,num2str(i),', ',num2str(ri),...
        ')']};
    title(str2);
    xlabel ('rp ')
    ylabel ('f r p aux')
%     disp (frp_aux')
    pause;
    close;
end
end % End function plot_2nd_int

function plot_3th_int (r, fr_aux)
% Descrpition: plot 3th integral

plot (r, fr_aux);
% str3 = {['3ª INT. (pos, valor): R desej max (eixo horizontal) = (',num2str(length(r)),', ',num2str(r(length(r))),')']};
str3 = {'3rd Integral'};
title(str3);
xlabel ('r ')
ylabel ('fr aux')
pause;
close;

end  % End function plot_3th_int