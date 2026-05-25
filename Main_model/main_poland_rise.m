%==========================================================================
% main_poland_rise.m
% --------------------------------------------------------------------------
% Full Bayesian estimation of NK MS-DSGE model for Poland using RISE.

clear; close all; clc;

cd("RISE_toolbox-master/")
rise_startup
cd ..

%==========================================================================
%% BLOCK 1: DATA PREPARATION
%==========================================================================
data_raw  = readtable('Data/Data_Poland.xlsx');

real_gdp  = data_raw.real_gdp;
deflator  = data_raw.gdp_deflator;
R_raw     = data_raw.Interest_rate;
b_raw     = data_raw.Debt_GDP;
s_raw     = data_raw.PrimSurplus;
dates_raw = data_raw.dates;
covid_raw = data_raw.covid;

pi_data    = log(deflator(2:end)) - log(deflator(1:end-1));
R_data = R_raw(2:end)/400;  
b_data = b_raw(2:end)/100;          
s_data = s_raw(2:end)/400;

log_gdp = log(real_gdp);
[trend, y_data_hp]    = hpfilter(log_gdp);

pi_target = 0.025 / 4;
R_target  =  mean(R_data);
s_target  = mean(s_data);

b_target = 0.5;

y_obs = y_data_hp(2:end);
pi_obs = pi_data  - pi_target;
R_obs  = R_data   - R_target;
b_obs  = b_data - b_target;
g_obs = -s_data + s_target;
covid_obs = covid_raw(2:end);

yt_mat = [y_obs, pi_obs, R_obs, g_obs, b_obs, covid_obs];

keep   = all(~isnan(yt_mat), 2);
yt_mat = yt_mat(keep, :);

start_date = '2000Q2';

db = ts(start_date, yt_mat, {'y','pi','R','g','b', 'dummy_covid'});

m=rise('nk_poland.rs');

startDate = datetime(2000,4,1);
endDate   = datetime(2025,10,1);

dates = (startDate:calmonths(3):endDate)';

%==========================================================================
%% BLOCK 2: DECLARING PRIORS
%==========================================================================

m = set(m, 'parameters', {
    'sigma',           1;
    'betta',           0.99;
    'kappa',           0.05;
    'delta_b_snk_2',   0;
    'pi_ss',           0.00625;
    'dy',              0.0089;
    'b_ss',            0.5;
});

priors = struct();

priors.rho_r      = {0.75, 0.75, 0.05, 'beta', 0.3, 0.98};
priors.rho_u      = {0.6, 0.6, 0.08, 'beta', 0, 0.95};
priors.rho_g      = {0.60, 0.60, 0.10, 'beta', 0.0, 0.98};
priors.rho_zeta   = {0.80, 0.80, 0.08, 'beta', 0.3, 0.99};
%priors.b_ss       = {0.4, 0.4, 0.05, 'normal'};
priors.R_ss       = {0.010,0.010,0.002,'normal'};
priors.g_ss       = {0.005,0.005,0.003,'normal'};
priors.gamma      = {0.08, 0.08, 0.03, 'normal', 0, 0.5};
priors.delta_y    = {0.05, 0.05, 0.03, 'normal', 0, 0.5};

priors.phi_pi_snk_1 = {1.8, 1.8, 0.1, 'normal', 1.1, 3.0};
priors.phi_pi_snk_2 = {0.8, 0.8, 0.1,'normal',0.5,1.0};
priors.phi_y_snk_1 = {0.30, 0.30, 0.05, 'normal', 0, 0.6};
priors.phi_y_snk_2  = {0.10,0.10,0.05,'normal',0,0.3};
priors.delta_b_snk_1 = {0.08,0.08,0.03,'gamma',0.01,0.5};

priors.sig_r        = {0.01,0.01,0.005,'inv_gamma'};
priors.sig_zeta     = {0.01,0.01,0.005,'inv_gamma'};
priors.sig_g        = {0.01, 0.01, 0.005, 'inv_gamma'};
priors.sig_u        = {0.01, 0.01, 0.003, 'inv_gamma'};
priors.sig_b        = {0.01,0.01,0.005,'inv_gamma'};
priors.sig_covid_y  = {0.03,0.03,0.01,'inv_gamma'};
priors.sig_covid_pi = {0.02,0.02,0.01,'inv_gamma'};

priors.snk_tp_1_2 = {0.04,0.04,0.015,'beta'};
priors.snk_tp_2_1 = {0.13,0.13,0.03,'beta'};

m=set(m,'data',db);
mest = estimate(m,'estim_priors',priors);

%%
smoothed_mest = filter(mest);

P_mest = smoothed_mest.smoothed_regime_probabilities;

p1_mest = P_mest.regime_1;
p2_mest = 1 - p1_mest;

Y_mest = [p1_mest p2_mest];

figure;
h = area(dates,Y_mest);

h(1).FaceColor = [0.75 0.75 0.75];
h(2).FaceColor = [0.35 0.35 0.35];

ylim([0 1])
ylabel('Probability')
grid on
%==========================================================================
%% BLOCK 3: MCMC WITH 200000 DRAWS
%==========================================================================
[objective,lb,ub,x0,SIG]=pull_objective(mest);
scale=0.0794*6;
myOpts=struct();
myOpts.tunedCov=scale*SIG;
myOpts.N=100000;
energy=@(varargin)-objective(varargin{:});
results=sample(rsamplers.rwmh(energy,x0,lb,ub,myOpts));


mddobj=mdd(results,energy,lb,ub,[],[],true);
bridge(mddobj,true,mdd.global_options);
%==========================================================================
%% BLOCK 4: EXTRACTING POSTERIOR DATA
%==========================================================================

burn = myOpts.N/10;
draws_struct = results{1,1}.pop;
draws_struct = draws_struct(burn+1:end);

N = length(draws_struct);
npar = length(draws_struct(1).x);

Theta = zeros(N,npar);

for n = 1:N
    Theta(n,:) = draws_struct(n).x(:)';
end

theta_mean = mean(Theta,1);
theta_median = median(Theta,1);

theta_lb90 = quantile(Theta,0.05,1);
theta_ub90 = quantile(Theta,0.95,1);

theta_lb68 = quantile(Theta,0.16,1);
theta_ub68 = quantile(Theta,0.84,1);

%%

param_names = { ...
    'rho_r'
    'rho_u'
    'rho_g'
    'rho_zeta'
    'b_ss'
    'R_ss'
    'g_ss'
    'gamma'
    'delta_y'
    'phi_pi_snk_1'
    'phi_pi_snk_2'
    'phi_y_snk_1'
    'phi_y_snk_2'
    'delta_b_snk_1'
    'sig_r'
    'sig_zeta'
    'sig_u'
    'sig_g'
    'sig_b'
    'sig_covid_y'
    'sig_covid_pi'
    'snk_tp_1_2'
    'snk_tp_2_1'
};

params_mean = struct();
params_ub = struct();
params_lb = struct();

T = table( ...
    theta_mean(:), ...
    theta_median(:), ...
    theta_lb68(:), ...
    theta_ub68(:), ...
    theta_lb90(:), ...
    theta_ub90(:), ...
    'RowNames', param_names, ...
    'VariableNames', { ...
        'Mean','Median', ...
        'LB68','UB68', ...
        'LB90','UB90'});

for i = 1:length(param_names)
    params_mean.(param_names{i}) = theta_mean(i);
end
mmean = set(m,'parameters',params_mean);

%==========================================================================
%% BLOCK 5: SMOOTHED POSTERIOR PROBABILITIES
%==========================================================================
smoothed_mean = filter(mmean);

P_mean = smoothed_mean.smoothed_regime_probabilities;

p1_mean = P_mean.regime_1;
p2_mean = 1 - p1_mean;

Y_mean = [p1_mean p2_mean];

figure;
h = area(dates,Y_mean);

h(1).FaceColor = [0.75 0.75 0.75];
h(2).FaceColor = [0.35 0.35 0.35];

ylim([0 1])
ylabel('Probability')
grid on

%==========================================================================
%% BLOCK 6: BAYESIAN REGIME-SPECIFIC IRFS
%==========================================================================

H = 40;

ndraws_irf = N;

idx = round(linspace(1,N,ndraws_irf));

shock_list = {'eps_zeta','eps_u','eps_MP','eps_g'};

var_list = {'y','pi','R','g','b'};

nshocks = length(shock_list);
nvars   = length(var_list);

nregimes = 2;
IRFS = zeros( ...
    ndraws_irf, ...
    H+1, ...
    nvars, ...
    nshocks, ...
    nregimes);

for j = 1:ndraws_irf

    theta_j = Theta(idx(j),:);

    params_j = struct();

    for i = 1:length(param_names)

        params_j.(param_names{i}) = theta_j(i);

    end

    m_j = set(m,'parameters',params_j);

    irfs_j = irf( ...
        m_j, ...
        'irf_shock_list',shock_list, ...
        'irf_var_list',var_list, ...
        'irf_periods',H);

    for s = 1:nshocks

        shockname = shock_list{s};

        for v = 1:nvars

            varname = var_list{v};

            tmp = irfs_j.(shockname).(varname);

            for r = 1:nregimes

                IRFS(j,:,v,s,r) = tmp(:,r)';

            end
        end
    end
end

IRF_MED  = median(IRFS,1);

IRF_LB68 = quantile(IRFS,0.16,1);
IRF_UB68 = quantile(IRFS,0.84,1);

IRF_LB90 = quantile(IRFS,0.05,1);
IRF_UB90 = quantile(IRFS,0.95,1);

%==========================================================================
%% PLOT REGIME-SPECIFIC IRFS
%==========================================================================

x = 0:H;

for r = 1:nregimes

    figure( ...
        'Name',['Regime ',num2str(r)], ...
        'Position',[100 100 1800 1000]);

    tiledlayout(nshocks,nvars, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for s = 1:nshocks

        shockname = shock_list{s};

        for v = 1:nvars

            varname = var_list{v};

            nexttile

            med = squeeze( ...
                IRF_MED(1,:,v,s,r));

            lb68 = squeeze( ...
                IRF_LB68(1,:,v,s,r));

            ub68 = squeeze( ...
                IRF_UB68(1,:,v,s,r));

            lb90 = squeeze( ...
                IRF_LB90(1,:,v,s,r));

            ub90 = squeeze( ...
                IRF_UB90(1,:,v,s,r));

            fill([x fliplr(x)], ...
                 [lb90 fliplr(ub90)], ...
                 [0.85 0.85 0.85], ...
                 'EdgeColor','none');

            hold on

            fill([x fliplr(x)], ...
                 [lb68 fliplr(ub68)], ...
                 [0.65 0.65 0.65], ...
                 'EdgeColor','none');

            plot(x,med,'k','LineWidth',2)

            yline(0,'k--')

            grid on
            if s == 1
                title(varname,'FontWeight','bold')
            end

            if v == 1

                ylabel( ...
                    ['Reaction to ', ...
                    strrep(shockname,'eps_',''), ...
                    ' shock'], ...
                    'FontWeight','bold');

            end

        end
    end

    sgtitle(['Posterior IRFs - Regime ',num2str(r)], ...
        'FontSize',16, ...
        'FontWeight','bold');

end


%==========================================================================
%% BLOCK 7: SMOOTHED POSTERIOR PROBABILITIES
%==========================================================================

smoothed_mean = filter(mmean);

P_mean = smoothed_mean.smoothed_regime_probabilities;

p1_mean = P_mean.regime_1;
p2_mean = 1 - p1_mean;

Y_mean = [p1_mean p2_mean];

figure;
h = area(dates,Y_mean);

h(1).FaceColor = [0.75 0.75 0.75];
h(2).FaceColor = [0.35 0.35 0.35];

ylim([0 1])
ylabel('Probability')
grid on

%==========================================================================
%% BLOCK 8: MODEL FIT VS DATA
%==========================================================================

hd_mean = historical_decomposition(mmean);
vars = {'y','pi','R','g','b'};


data_series = struct();

data_series.y  = y_obs;
data_series.pi = pi_obs;
data_series.R  = R_obs;
data_series.g  = g_obs;
data_series.b  = b_obs;

model_series = struct();

for i = 1:length(vars)

    varname = vars{i};

    hd_obj = hd_mean.(varname);

    init_col = strcmp(hd_obj.varnames,'init');

    shock_part = sum( ...
        hd_obj.data(:,~init_col), ...
        2);

    init_part = hd_obj.data(:,init_col);

    model_series.(varname) = ...
        shock_part + init_part;

end
model_series.b = model_series.b + 0.07;

dates_fit = dates;

figure( ...
    'Name','Model Fit vs Data', ...
    'Position',[100 100 1400 800]);

tiledlayout(3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

nexttile

plot(dates_fit, ...
    data_series.y, ...
    'k', ...
    'LineWidth',2);

hold on

plot(dates_fit, ...
    model_series.y, ...
    '--r', ...
    'LineWidth',2);

title('Output Gap')

legend('Data','Model', ...
    'Location','best')

grid on

nexttile

plot(dates_fit, ...
    data_series.pi, ...
    'k', ...
    'LineWidth',2);

hold on

plot(dates_fit, ...
    model_series.pi, ...
    '--r', ...
    'LineWidth',2);

title('Inflation')

legend('Data','Model', ...
    'Location','best')

grid on

nexttile

plot(dates_fit, ...
    data_series.R, ...
    'k', ...
    'LineWidth',2);

hold on

plot(dates_fit, ...
    model_series.R, ...
    '--r', ...
    'LineWidth',2);

title('Interest Rate')

legend('Data','Model', ...
    'Location','best')

grid on
nexttile

plot(dates_fit, ...
    data_series.g, ...
    'k', ...
    'LineWidth',2);

hold on

plot(dates_fit, ...
    model_series.g, ...
    '--r', ...
    'LineWidth',2);

title('Government Spending')

legend('Data','Model', ...
    'Location','best')

grid on

nexttile

plot(dates_fit, ...
    data_series.b, ...
    'k', ...
    'LineWidth',2);

hold on

plot(dates_fit, ...
    model_series.b, ...
    '--r', ...
    'LineWidth',2);

title('Debt')

legend('Data','Model', ...
    'Location','best')

grid on

sgtitle( ...
    'Model Fit: Historical Decomposition Reconstruction vs Data', ...
    'FontSize',16, ...
    'FontWeight','bold');
