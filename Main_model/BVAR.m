clear; clc; close all; rng(1);

%% ==========================================================
% LOAD DATA
%% ==========================================================

[num,txt,raw] = xlsread('Data/Data_final.xlsx');

fiscal_ratio_raw = cell2mat(raw(2,2:end));
nominal_rate_raw = cell2mat(raw(3,2:end));
real_gdp_raw     = cell2mat(raw(4,2:end));
deflator_raw     = cell2mat(raw(5,2:end));

dates_raw = raw(1,2:end);

Time = length(dates_raw);
dates_dt = NaT(Time,1);

for i = 1:Time
    str = dates_raw{i};
    year = str2double(str(1:4));
    quarter = str2double(str(end));
    month = (quarter-1)*3 + 1;
    dates_dt(i) = datetime(year,month,1);
end

%% ==========================================================
% TRANSFORMATIONS
% Inflation, output growth = YoY log differences
%% ==========================================================

inflation  = 100*(log(deflator_raw(5:end)) - log(deflator_raw(1:end-4)));
gdp_growth = 100*(log(real_gdp_raw(5:end)) - log(real_gdp_raw(1:end-4)));

nominal_rate = nominal_rate_raw(5:end);
fiscal_ratio = fiscal_ratio_raw(5:end);

real_rate = nominal_rate - inflation;

dates_dt = dates_dt(5:end);

%% ==========================================================
% DATA MATRIX
% 1 inflation
% 2 output growth
% 3 real rate
% 4 fiscal ratio
%% ==========================================================

Z = [inflation' gdp_growth' real_rate' fiscal_ratio'];

muZ = mean(Z);
sdZ = std(Z);
Z = (Z - muZ)./sdZ;

[T,n] = size(Z);

Y = Z(2:end,:);
X = Z(1:end-1,:);
X = [ones(size(X,1),1) X];

T_eff = size(Y,1);
dimX  = size(X,2);

dates_var = dates_dt(2:end);

%% ==========================================================
% SETTINGS
%% ==========================================================

k = 2;
n_iter  = 100000;
burn_in = 10000;

%% ==========================================================
% MINNESOTA PRIOR
%% ==========================================================

lambda1 = 0.18;
lambda2 = 0.25;

sig2 = var(Z);

B0 = zeros(n,dimX);

for eq = 1:n
    B0(eq,1+eq) = 1;
end

V_prior = cell(n,1);

for eq = 1:n
    
    V = zeros(dimX);
    
    V(1,1) = 5^2;
    
    for j = 1:n
        
        col = 1 + j;
        
        if eq == j
            v = lambda1^2;
        else
            v = lambda1^2 * lambda2^2 * sig2(eq)/sig2(j);
        end
        
        V(col,col) = v;
    end
    
    V_prior{eq} = V;
end

%% ==========================================================
% REGIME-SPECIFIC INFORMATIVE PRIORS
%% ==========================================================

B0_r1 = B0;
B0_r2 = B0;

B0_r1(3,2) = 1.00;  
B0_r2(3,2) = 0.10;   

B0_r1(4,5) = 1.00;   
B0_r2(4,5) = 0.30;

%% ==========================================================
% PRIORS FOR SIGMA AND TRANSITIONS
%% ==========================================================

nu0 = n + 8;
S0  = 0.7*eye(n);

alpha0 = [300 1;
          3 100];
%% ==========================================================
% INITIAL VALUES
%% ==========================================================

Phi_fixed = zeros(2,dimX);

Phi_sw = cell(k,1);

Phi_sw{1} = B0_r1(3:4,:);
Phi_sw{2} = B0_r2(3:4,:);

Sigma = cell(k,1);
Sigma{1} = 0.6*eye(n);
Sigma{2} = 1.8*eye(n);

P = [0.98 0.02;
     0.02 0.98];

%% ==========================================================
% STORAGE
%% ==========================================================

P_store = zeros(k,k,n_iter);
prob_store = zeros(T_eff,n_iter);
Phi_store = zeros(n,dimX,k,n_iter);

%% ==========================================================
% GIBBS LOOP
%% ==========================================================

for iter = 1:n_iter
    
    %% ------------------------------------------------------
    % STEP 1: FIXED EQUATIONS (1,2)
    %% ------------------------------------------------------
    
    for eq = 1:2
        
        V0 = V_prior{eq};
        b0 = B0(eq,:)';
        
        Prec = inv(V0) + X'*X;
        Mean = Prec \ (inv(V0)*b0 + X'*Y(:,eq));
        VarP = inv(Prec);
        
        beta = Mean + chol(VarP,'lower')*randn(dimX,1);
        
        Phi_fixed(eq,:) = beta';
    end
    
    %% ------------------------------------------------------
    % STEP 2: HAMILTON FILTER
    %% ------------------------------------------------------
    
    xi_pred = zeros(T_eff,k);
    xi_filt = zeros(T_eff,k);
    like = zeros(T_eff,k);
    
    xi_pred(1,:) = [0.95 0.05];
    
    for t = 1:T_eff
        
        for j = 1:k
            
            Phi = zeros(n,dimX);
            Phi(1:2,:) = Phi_fixed;
            Phi(3:4,:) = Phi_sw{j};
            
            mu = Phi*X(t,:)';
            
            res = Y(t,3:4)' - mu(3:4);
            Sig = Sigma{j}(3:4,3:4);
            
            like(t,j) = (mvnpdf(res',zeros(1,2),Sig)+1e-12)^1.2;
        end
        
        if t==1
            prior = [0.95 0.05];
        else
            prior = xi_pred(t,:);
        end
        
        nume = prior .* like(t,:);
        xi_filt(t,:) = nume / sum(nume);
        
        if t<T_eff
            xi_pred(t+1,:) = xi_filt(t,:)*P;
        end
    end
    
    %% ------------------------------------------------------
    % STEP 3: SMOOTHER
    %% ------------------------------------------------------
    
    xi = zeros(T_eff,k);
    xi(end,:) = xi_filt(end,:);
    
    for t = T_eff-1:-1:1
        
        for i = 1:k
            
            temp = 0;
            
            for j = 1:k
                temp = temp + P(i,j)*xi(t+1,j) ...
                    / max(xi_pred(t+1,j),1e-12);
            end
            
            xi(t,i) = xi_filt(t,i)*temp;
        end
        
        xi(t,:) = xi(t,:)/sum(xi(t,:));
    end
    
    %% ------------------------------------------------------
    % STEP 4: DRAW SWITCHING EQUATIONS
    %% ------------------------------------------------------
    
    for j = 1:k
        
        W = diag(xi(:,j));
        
        if j==1
            PriorMean = B0_r1;
        else
            PriorMean = B0_r2;
        end
        
        eq = 3;
        V0 = V_prior{eq};
        b0 = PriorMean(eq,:)';
        
        Prec = inv(V0) + X'*W*X;
        Mean = Prec \ (inv(V0)*b0 + X'*W*Y(:,eq));
        VarP = inv(Prec);
        
        beta = Mean + chol(VarP,'lower')*randn(dimX,1);
        Phi_sw{j}(1,:) = beta';
        
        eq = 4;
        V0 = V_prior{eq};
        b0 = PriorMean(eq,:)';
        
        Prec = inv(V0) + X'*W*X;
        Mean = Prec \ (inv(V0)*b0 + X'*W*Y(:,eq));
        VarP = inv(Prec);
        
        beta = Mean + chol(VarP,'lower')*randn(dimX,1);
        Phi_sw{j}(2,:) = beta';
        
        %% Draw Sigma
        
        Phi = zeros(n,dimX);
        Phi(1:2,:) = Phi_fixed;
        Phi(3:4,:) = Phi_sw{j};
        
        resid = Y - X*Phi';
        
        S_post = S0;
        
        for t = 1:T_eff
            S_post = S_post + xi(t,j)*(resid(t,:)'*resid(t,:));
        end
        
        Sigma{j} = iwishrnd(S_post,nu0+sum(xi(:,j)));
    end
    
    %% ------------------------------------------------------
    % STEP 5: TRANSITION MATRIX
    %% ------------------------------------------------------
    
    for i = 1:k
        
        stay = sum(xi(1:end-1,i).*xi(2:end,i));
        move = sum(xi(1:end-1,i).*xi(2:end,3-i));
        
        P(i,:) = dirichlet_rnd(alpha0(i,:) + [stay move]);
        P_store(:,:,iter) = P;
    end
    
    %% ------------------------------------------------------
    % IDENTIFICATION:
    %% ------------------------------------------------------
    
    if Phi_sw{1}(1,2) < Phi_sw{2}(1,2)
        
        Phi_sw = Phi_sw([2 1]);
        Sigma  = Sigma([2 1]);
        P      = P([2 1],[2 1]);
        xi     = xi(:,[2 1]);
    end
    
    %% ------------------------------------------------------
    % STORE
    %% ------------------------------------------------------
    
    prob_store(:,iter) = xi(:,1);
    
    for j = 1:k
        
        Phi = zeros(n,dimX);
        Phi(1:2,:) = Phi_fixed;
        Phi(3:4,:) = Phi_sw{j};
        
        Phi_store(:,:,j,iter) = Phi;
    end
    
end

%% ==========================================================
% POSTERIOR RESULTS
%% ==========================================================

prob_regime1 = mean(prob_store(:,burn_in:end),2);

Phi1 = mean(Phi_store(:,:,1,burn_in:end),4);
Phi2 = mean(Phi_store(:,:,2,burn_in:end),4);

%% ==========================================================
% PLOT
%% ==========================================================

figure; hold on;

x = dates_var;
p = movmean(prob_regime1,4);

Yplot = [p , 1-p];

h = area(x,Yplot);

h(1).FaceColor = [0.75 0.75 0.75];
h(1).EdgeColor = 'none';

h(2).FaceColor = [0.35 0.35 0.35];
h(2).EdgeColor = 'none';

ylim([0 1]);
grid on;

ylabel('Probability');

xticks(x(1:4:end));
xtickformat('yyyy');

%% ==========================================================
% DISPLAY POLICY COEFFICIENTS
%% ==========================================================

disp('===== REGIME 1 (Monetary Dominance) =====')
disp(Phi1)

disp('===== REGIME 2 (Fiscal Dominance) =====')
disp(Phi2)

%% ==========================================================
% DISPLAY AVARAGES ACROSS REGIMES
%% ==========================================================

p1 = prob_regime1;
p2 = 1 - p1;

avg_regime1 = zeros(1,4);
avg_regime2 = zeros(1,4);

Zraw = [inflation' gdp_growth' real_rate' fiscal_ratio'];
Zraw = Zraw(2:end,:);

for j = 1:4
    avg_regime1(j) = sum(p1 .* Zraw(:,j)) / sum(p1);
    avg_regime2(j) = sum(p2 .* Zraw(:,j)) / sum(p2);
end

disp('===== REGIME 1 AVERAGES =====')
disp(avg_regime1)

disp('===== REGIME 2 AVERAGES =====')
disp(avg_regime2)

steady_states = zeros(2,n);

c1 = Phi1(:,1);       
A1 = Phi1(:,2:end);    
steady_states(1,:) = ((eye(n) - A1) \ c1)';

c2 = Phi2(:,1);       
A2 = Phi2(:,2:end);    
steady_states(2,:) = ((eye(n) - A2) \ c2)';

disp('===== STEADY STATES =====')
disp(steady_states)

P_post = mean(P_store(:,:,burn_in:end),3);

disp('===== POSTERIOR TRANSITION MATRIX =====')
disp(P_post)
%% ==========================================================
% FUNCTIONS
%% ==========================================================

function x = dirichlet_rnd(alpha)
y = gamrnd(alpha,1);
x = y/sum(y);
end