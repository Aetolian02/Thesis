%==========================================================================
% main_poland_rise.m
% --------------------------------------------------------------------------
% Full Bayesian estimation of NK MS-DSGE model for Poland using RISE.


clear; close all; clc;

fprintf('=============================================\n')
fprintf(' NK MS-DSGE Poland — RISE estimation\n')
fprintf('=============================================\n\n')
cd("RISE_toolbox-master/")
rise_startup
cd ..

%==========================================================================
%% BLOCK 1: DATA PREPARATION
%==========================================================================
fprintf('Block 1: Loading data...\n')

data_raw  = readtable('Data/Data_Poland.xlsx');

real_gdp  = data_raw.real_gdp;
deflator  = data_raw.gdp_deflator;
R_raw     = data_raw.Interest_rate;
b_raw     = data_raw.Debt_GDP;
s_raw     = data_raw.PrimSurplus;

pi_data    = log(deflator(2:end)) - log(deflator(1:end-1));
R_data = R_raw(2:end)/400;  
b_data = b_raw(2:end)/100;          
s_data = s_raw(2:end)/100;

% --- Output gap via HP filter ---
log_gdp = log(real_gdp);
[trend, y_data_hp]    = hpfilter(log_gdp);

% --- Output gap via detrending ---
y_data_lin = log(real_gdp(2:end)) - log(real_gdp(1:end-1));

y_target = mean(y_data_lin);
pi_target = 0.025 / 4;
R_target  = mean(R_data);
b_target = 0.55;
s_target  = mean(s_data);

%y_obs  = y_data_hp(2:end);
y_obs = y_data_hp(2:end);
pi_obs = pi_data  - pi_target;
R_obs  = R_data   - R_target;
b_obs  = b_data - b_target;
g_obs = -s_data + s_target;

yt_mat = [y_obs, pi_obs, R_obs, b_obs, g_obs];

keep   = all(~isnan(yt_mat), 2);
yt_mat = yt_mat(keep, :);
T      = sum(keep);

start_date = '2000Q2';

db = ts(start_date, yt_mat, {'y','pi','R','b','g'});

m=rise('nk_poland.rs');

startDate = datetime(2000,4,1);
endDate   = datetime(2025,10,1);

dates = (startDate:calmonths(3):endDate)';

%==========================================================================
%% BLOCK 2: DECLARING PRIORS
%==========================================================================
fprintf('Block 2: Building RISE model...\n')

m = set(m, 'parameters', {
    'sigma',           1.5;
    'betta',           0.99;
    'kappa',           0.05;
    'delta_b_snk_2',   0;
    'pi_ss',           0.0063;
    'dy',              0.01;
});

priors = struct();

priors.rho_r      = {0.7,  0.7,  0.05, 'beta',   0,    0.99};  
priors.rho_u      = {0.5,  0.5,  0.15, 'beta',   0,    0.99};
priors.rho_g      = {0.5,  0.5,  0.15, 'beta',   0,    0.99};
priors.rho_zeta   = {0.8,  0.8,  0.15, 'beta',   0,    0.99};
priors.b_ss       = {0.55, 0.55, 0.03, 'normal'};               
priors.R_ss       = {0.01, 0.01, 0.003,'normal'};              
priors.g_ss       = {0.02, 0.02, 0.005,'normal'};
priors.gamma      = {0.1,  0.1,  0.005,'normal', 0,    0.99};

priors.psi_pi_snk_1  = {1.8, 1.8, 0.20, 'normal', 1.01, 3.0};  
priors.psi_pi_snk_2  = {0.7, 0.7, 0.15, 'normal', 0.20, 0.99};
priors.delta_b_snk_1 = {0.05,0.05,0.02, 'gamma',  0.01, 0.50};

priors.sig_r    = {0.02, 0.02, 0.01, 'inv_gamma'};
priors.sig_zeta = {0.02, 0.02, 0.01, 'inv_gamma'};
priors.sig_u    = {0.02, 0.02, 0.01, 'inv_gamma'};
priors.sig_g    = {0.02, 0.02, 0.01, 'inv_gamma'};
priors.sig_b    = {0.02, 0.02, 0.01, 'inv_gamma'};

priors.snk_tp_1_2 = {0.02, 0.02, 0.01, 'beta'};
priors.snk_tp_2_1 = {0.12, 0.12, 0.03, 'beta'};

m=set(m,'data',db);
mest = estimate(m,'estim_priors',priors);

%==========================================================================
%% BLOCK 3: CREATING IRFS ACROSS REGIMES
%==========================================================================
mest_irfs = irf(mest,'irf_periods',400);
quick_irfs(mest,mest_irfs,get(mest,'endo_list(original)'))

%==========================================================================
%% BLOCK 4: SMOOTHED POSTERIOR PROBABILITIES
%==========================================================================
smoothed = filter(mest);

P = smoothed.smoothed_regime_probabilities;

p1 = P.regime_1;
p2 = 1 - p1;

Y = [p1 p2];

figure;
h = area(dates,Y);

h(1).FaceColor = [0.75 0.75 0.75];
h(2).FaceColor = [0.35 0.35 0.35];

ylim([0 1])
ylabel('Probability')
title('Regime Probabilities')
grid on


%==========================================================================
%% BLOCK 5: HISTORICAL DECOMPOSITION: TO BE DONE
%==========================================================================
hd = historical_decomposition(mest);

%==========================================================================
%% BLOCK 6: MCMC WITH 100000 DRAWS
%==========================================================================
[objective,lb,ub,x0,SIG]=pull_objective(mest);
scale=0.0794;
myOpts=struct();
myOpts.tunedCov=scale*SIG;
myOpts.N=100000;
energy=@(varargin)-objective(varargin{:});
results=sample(rsamplers.rwmh(energy,x0,lb,ub,myOpts));


mddobj=mdd(results,energy,lb,ub,...
[],[],true);
bridge(mddobj,true,mdd.global_options)