@endogenous
    y, pi, R, g, b, zeta, u, dummy_covid
;

@exogenous
    eps_zeta, eps_u, eps_MP, eps_g, eps_b, eps_covid
;

@parameters
    sigma, betta, kappa, rho_r, rho_zeta, rho_u, b_ss, rho_g, gamma, dy, R_ss, pi_ss, delta_y, g_ss,
    sig_zeta, sig_u, sig_g, sig_r, sig_b, sig_covid_y, sig_covid_pi,
    snk_tp_1_2, snk_tp_2_1
;

@parameters(snk,2)
    delta_b, phi_pi, phi_y
;
@observables 
    y, pi, R, g, b, dummy_covid
;

@model
    y{t} = y{t+1} - (1/sigma)*(R{t} - pi{t+1}) + gamma*(g{t}-g{t+1}) + zeta{t} - sig_covid_y*dummy_covid{t};
    pi{t} = betta*(1+dy)^(sigma)*pi{t+1} + kappa*y{t} + u{t} + sig_covid_pi*dummy_covid{t};
    R{t} = rho_r*R{t-1} + (1-rho_r)*(phi_pi*pi{t}+phi_y*y{t}) + sig_r*eps_MP{t};
    b{t} = (1+R_ss-pi_ss-dy)*b{t-1} + (R_ss+R{t-1}-pi{t}-pi_ss-dy-y{t}+y{t-1})*b_ss + g{t} + g_ss + sig_b*eps_b{t};
    zeta{t} = rho_zeta*zeta{t-1}+sig_zeta*eps_zeta{t};
    u{t} = rho_u*u{t-1}+sig_u*eps_u{t};
    g{t} = rho_g*g{t-1}-(1-rho_g)*(delta_b*b{t-1} + delta_y*y{t}) + sig_g*eps_g{t};
    dummy_covid{t} = eps_covid{t};
