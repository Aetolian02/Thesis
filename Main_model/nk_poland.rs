@endogenous
    y, pi, R, b, zeta, u, g
;

@exogenous
    eps_zeta, eps_u, eps_MP, eps_g, eps_b
;

@parameters
    sigma, betta, kappa, rho_r, rho_zeta, rho_u, b_ss, rho_g, gamma, dy, R_ss, pi_ss, g_ss
    sig_zeta, sig_u, sig_g, sig_r, sig_b,
    snk_tp_1_2, snk_tp_2_1
;

@parameters(snk,2)
    delta_b, psi_pi
;
@observables 
    y, pi, R, b, g
;

@model
    y{t} = y{t+1} - (1/sigma)*(R{t} - pi{t+1}) + gamma*(g{t}-g{t+1}) + zeta{t};
    pi{t} = betta*(1+dy)^(1/sigma)*pi{t+1} + kappa*y{t} + u{t};
    R{t} = rho_r*R{t-1} + (1-rho_r)*psi_pi*pi{t} + sig_r*eps_MP{t};
    b{t} = (1+R_ss+R{t-1}-pi{t}-pi_ss-dy-y{t}+y{t-1})*b{t-1} + (R_ss+R{t-1}-pi{t}-pi_ss-dy-y{t}+y{t-1})*b_ss + g{t} + g_ss + sig_b*eps_b;
    zeta{t} = rho_zeta*zeta{t-1}+sig_zeta*eps_zeta{t};
    u{t} = rho_u*u{t-1}+sig_u*eps_u{t};
    g{t} = rho_g*g{t-1}-delta_b*b{t-1}+sig_g*eps_g{t};