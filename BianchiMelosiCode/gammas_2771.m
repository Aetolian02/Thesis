%==================================================================
% Replication codes for
% The Dire Effects of the Lack of Monetary and Fiscal Coordination, 
% by Francesco Bianchi and Leonardo Melosi
% Journal of Monetary Economics, Volume 104, pp. 1-22. 
%==================================================================


function [GAM0,GAM1,PSI,PPI,snames] = gammas_2771(p_st,p_ss,msel,varargin)
redD=msel(22);
nst=msel(10);
if nst>1
    if nargin>3 
        HB=varargin{1};
        ZLB=0;
    else 
        ZLB=1; 
    end
end

%---------------------------------
psi_zlb    = p_st(1,1);
kap        = p_st(2,1);
psipi      = p_st(3,1);
psiy       = p_st(4,1);
rhoR       = p_st(5,1);
rhochi     = p_st(6,1);if rhochi~=0;aaaaaaa;end
scaleMPsh  = p_st(7,1);
deltab     = p_st(8,1);
deltae     = p_st(9,1);
rhotau     = p_st(10,1);

% rhodem     = p_st(11,1);
rhoe_ST    = p_st(12,1);
psiDy      = p_st(13,1);
deltay0    = p_st(14,1);
deltay1    = p_st(15,1);
varsigma   = p_st(16,1); % 
Cap_phi    = p_st(17,1); % 
varphi     = p_st(18,1); 
deltaestar = p_st(19,1);

% Autocorrelations for the shocks                
rhozetaMP  = p_st(21,1);                
rhog       = p_st(22,1);   
rhozetaA   = p_st(23,1);    
rhozetaT   = p_st(24,1);
rhozetaD   = p_st(25,1);
rhozetaES  = p_st(26,1);
rhozetaTP  = p_st(27,1);
rhozetaMU  = p_st(28,1);
rhozetaEL  = p_st(29,1);
rhozetaXX  = p_st(30,1);


% Intercepts for the shocks                
czetaMP  = p_st(31,:);               
czetag   = p_st(32,:);     if czetag~=0;ddddddd;end
czetaA   = p_st(33,:);     if czetaA~=0;ddddddd;end
czetaT   = p_st(34,:);     if czetaT~=0;ddddddd;end
czetaD   = p_st(35,:);
czetaES  = p_st(36,:);     if czetaES~=0;ddddddd;end
czetaTP  = p_st(37,:);     if czetaTP~=0;ddddddd;end
czetaMU  = p_st(38,:);     if czetaMU~=0;ddddddd;end
czetaLT  = p_st(39,:);     if czetaLT~=0;ddddddd;end
czetaXX  = p_st(40,:);     if czetaXX~=0;ddddddd;end

if nst>1 && ZLB==1
    if sum(p_st(41:50,:)~=0)~=sum(p_st(31:40,:)~=0);ddddddddd;end              
    EczetaMP  = p_st(41,1);                
    Eczetachi = p_st(42,1);   
    EczetaA   = p_st(43,1);    
    EczetaT   = p_st(44,1);
    EczetaD   = p_st(45,1);
    EczetaES  = p_st(46,1);
    EczetaTP  = p_st(47,1);
    EczetaMU  = p_st(48,1);
    EczetaLT  = p_st(49,1);
    EczetaXX  = p_st(50,1);
end


%-----------------------------------
% rr_st      = p_ss(1);
pi_st       = p_ss(2);
gam_st      = p_ss(3);
b_st        = p_ss(4);
g_st        = p_ss(5);
tau_st      = p_ss(6);

iota_y      = p_ss(10);if iota_y~=0;g is assumed exogenous;end
phi_y       = p_ss(11);

beta        = p_ss(14);
ave_mat     = p_ss(15);
alpha       = p_ss(16);
iota_e      = p_ss(17);if iota_e~=0;g is assumed exogenous;end

rr_steady = (1+gam_st) / beta-1;

R_st     = 1+(pi_st + rr_steady);
beta1 = 1/beta;
rhowood=(1-1/ave_mat)*beta1;

if rhowood>beta1;aaddda; end

%------------------------------------------------------------------------
snames=cell(1,200);
ii=1;
vy    = ii;ii=ii+1;snames{ii-1}=('y');%  [1] Output
vpi   = ii;ii=ii+1;snames{ii-1}=('\pi');%  [2] Inflation
vR    = ii;ii=ii+1;snames{ii-1}=('R');%  [3] FFR
vb    = ii;ii=ii+1;snames{ii-1}=('b');%  [4] Debt

vg    = ii;ii=ii+1;snames{ii-1}=('g');%  [5] g (transormation of G/Y)
vtau  = ii;ii=ii+1;snames{ii-1}=('\tau'); %[6]%  [6] Taxes
va    = ii;ii=ii+1;snames{ii-1}=('z');%  [7] TFP
ve    = ii;ii=ii+1;snames{ii-1}=('tr');%  [8]Total transfers

vmupR = ii;ii=ii+1;snames{ii-1}=('\mu');% [9] Mark-up shock (rescaled)
vtp   = ii;ii=ii+1;snames{ii-1}=('tp');%   [10]Term premium 
ve_LT = ii;ii=ii+1;snames{ii-1}=('e^L');%   [11]Long term comp of transfers
vynat = ii;ii=ii+1;snames{ii-1}=('y^*');%     [12] Natural output

vdem  = ii;ii=ii+1;snames{ii-1}=('d');%      [13] Demand shock
vQm   = ii;ii=ii+1;snames{ii-1}=('Q');%     [14] Price long term bond
vRLT  = ii;ii=ii+1;snames{ii-1}=('R^L');%     [15] Return long term bond

% Expectation variables
jj=1;
vEy   = ii;ii=ii+1;jj=jj+1;snames{ii-1}=('E(y)'); %   [16] Expected output
vEpi  = ii;ii=ii+1;jj=jj+1;snames{ii-1}=('E(\pi)');%    [17] Expected inflation
vERLT = ii;ii=ii+1;jj=jj+1;snames{ii-1}=('E(R^L)');%    [19] Expected return maturity bond 

%Dummy Variables
if nst>1
    if ZLB==0; % All dummy variables inside gamma file
        ves1       = ii;snames{ii-1}='ves1';ii=ii+1;
        ii=ii+nst-1;
        vesEND=ves1+nst-1;
    else % Dummy variables modelled using MSDDGE setup
        ves1       = ii;snames{ii-1}='ves1';ii=ii+1;
        vesEND=ves1;
    end
end

%------------------------------
nend = jj-1;
neq  = ii-1;

snames=snames(:,1:neq);

%------------------------------
% Exog shocks
ii=1;
eR     = ii;ii=ii+1; % [1] MP
eg     = ii;ii=ii+1; % [2] g
ea     = ii;ii=ii+1; % [3] tfp

etau   = ii;ii=ii+1; % [4] taxation
edem   = ii;ii=ii+1; % [5] demand shock
ee_ST  = ii;ii=ii+1; % [6] short term component tr

etp    = ii;ii=ii+1; % [7] term premium
e_mupR = ii;ii=ii+1; % [8] mark-up
ee_LT  = ii;ii=ii+1; % [9] long term component tr

if nst>1
    sh_es1       = ii;   ii=ii+1;
    if ZLB==0
        ii=ii+nst-1; 
    end
end

nex  = ii-1;

%------------------------------
GAM0 = zeros(neq,neq);
GAM1 = zeros(neq,neq);
PSI  = zeros(neq,nex);
PPI  = zeros(neq,nend);

%-----------------------------
% Define some composite parameters
Cap_phi_g1 = Cap_phi/(1+gam_st);
eq=1;

%-----------------------------
% [1] IS Curve
GAM0(eq,vy)   = 1+Cap_phi_g1;
GAM0(eq,vR)   = 1-Cap_phi_g1;
GAM0(eq,vg)   = -(1+Cap_phi_g1-rhog);
GAM0(eq,va)   = -(rhozetaA-Cap_phi_g1);
GAM0(eq,vEy)  = -1;
GAM0(eq,vEpi) = -(1-Cap_phi_g1);
if redD==2 % Structural form for preference shock
    GAM0(eq,vdem) = -(1-Cap_phi_g1)*(1-rhozetaD);
elseif redD==1 || redD==-1% Reduced form for preference shock
    GAM0(eq,vdem) = -1;
end

GAM1(eq,vy)   = Cap_phi_g1;
GAM1(eq,vg)   = -Cap_phi_g1;
if nst>1
    if redD~=-1 % Discrete demand shock is separate from continuous demand shock
        coeffC=-(1-Cap_phi_g1)*czetaD;
        if ZLB==0
            coeffE=((1-Cap_phi_g1)*czetaD)*HB;
        else
            coeffE=(+(1-Cap_phi_g1)*EczetaD);
        end
        GAM0(eq,ves1:vesEND)    = coeffE+coeffC;
    end
end
eq=eq+1;


%-----------------------------
% [2] Phillips Curve   
varsigbeta_p_1 = 1+varsigma*beta;
GAM0(eq,vy)    = -(kap/varsigbeta_p_1)*(  1/(1-Cap_phi_g1)+(varphi+alpha)/(1-alpha) );
GAM0(eq,vpi)   = 1;
GAM0(eq,vmupR) = -1;
GAM0(eq,vg)    = kap/((1-Cap_phi_g1)*varsigbeta_p_1);
GAM0(eq,va)    = -kap*Cap_phi_g1/((1-Cap_phi_g1)*varsigbeta_p_1);

GAM0(eq,vEpi)  = -beta/varsigbeta_p_1;

GAM1(eq,vy)    = -kap*Cap_phi_g1/((1-Cap_phi_g1)*varsigbeta_p_1);
GAM1(eq,vpi)   = varsigma/varsigbeta_p_1;
GAM1(eq,vg)    = kap*Cap_phi_g1/((1-Cap_phi_g1)*varsigbeta_p_1); 
eq=eq+1;


%-----------------------------
% [3] Monetary policy rule 
GAM0(eq,vR)      = 1;
GAM0(eq,vy)      = -(psiy+psiDy)*(1-rhoR);
GAM0(eq,va)      = -psiDy*(1-rhoR);
GAM0(eq,vynat)   = psiy*(1-rhoR);
GAM0(eq,vpi)     = -psipi*(1-rhoR);
if nst>1;GAM0(eq,ves1:vesEND) = (1-rhoR)*psi_zlb*(R_st-1);end


GAM1(eq,vR)      = rhoR;
GAM1(eq,vy)      = -psiDy*(1-rhoR);

PSI(eq,eR)    = scaleMPsh;
eq=eq+1;

%-----------------------------
% [4] G
GAM0(eq,vg)    = 1;
GAM1(eq,vg)    = rhog;
PSI(eq,eg)     = 1;
eq=eq+1;

%-----------------------------
% [5] Tax Rule  
GAM0(eq,vtau)    = 1;
GAM0(eq,vy)      = -deltay0*(1-rhotau);
GAM0(eq,vynat)   = deltay0*(1-rhotau);

GAM1(eq,vtau)    = rhotau;
GAM1(eq,vy)      = deltay1*(1-rhotau);
GAM1(eq,vynat)   = -deltay1*(1-rhotau);
GAM1(eq,vb)      = deltab*(1-rhotau);

GAM0(eq,ve)      = -deltae*(1-rhotau);
GAM0(eq,ve_LT)   = -deltaestar*(1-rhotau);
GAM0(eq,vg)      = -(deltae+deltaestar)*(1-rhotau)/g_st;
PSI(eq,etau)     = 1;
eq=eq+1;

%-----------------------
% [6] Debt
bc1=b_st*beta1;
GAM0(eq,vy)     = bc1;
GAM0(eq,vpi)    = bc1;
GAM0(eq,vtau)   = 1;      
GAM0(eq,vb)     = 1;
GAM0(eq,va)     = bc1;
GAM0(eq,ve)     = -1;      
GAM0(eq,vRLT)   = -bc1;    
GAM0(eq,vg)     = -1/g_st; 
GAM0(eq,vtp)    = -1;      

GAM1(eq,vb)     = beta1;
GAM1(eq,vy)     = bc1;
eq=eq+1;

%-----------------------
% [7] Return Long term bond
GAM0(eq,vRLT)    = 1; 
GAM0(eq,vQm)     = -rhowood*(R_st^-1); 
GAM1(eq,vQm)     = -1; 
eq=eq+1;

%-----------------------
% [8] No arbitrage
GAM0(eq,vR)     = 1; 
GAM0(eq,vERLT)  = -1;
eq=eq+1;

%-----------------------
% [9] total transfers
GAM0(eq,ve)     = 1;
GAM0(eq,vy)     = -(1-rhoe_ST)*phi_y;
GAM0(eq,vynat)  = (1-rhoe_ST)*phi_y;
GAM0(eq,ve_LT)  = -(1-rhoe_ST);

GAM1(eq,ve)     = rhoe_ST;
PSI(eq,ee_ST)   = 1;
eq=eq+1;

%-----------------------
% [10] Long term component transfers
GAM0(eq,ve_LT)  = 1;
GAM1(eq,ve_LT)  = rhozetaEL;
PSI(eq,ee_LT)   = 1;
eq=eq+1;

%-----------------------
% [11] TP
GAM0(eq,vtp) = 1;
GAM1(eq,vtp) = rhozetaTP;
PSI(eq,etp)  = 1;
eq=eq+1;

%-----------------------------
% [12] Technology Growth
GAM0(eq,va)    = 1;
GAM1(eq,va)    = rhozetaA;
PSI(eq,ea)     = 1;
eq=eq+1;

%-----------------------
% [13] Demand shock
GAM0(eq,vdem) = 1;
if nst>1 && redD==-1;GAM0(eq,ves1:vesEND)  =  -czetaD;end
GAM1(eq,vdem) = rhozetaD;
PSI(eq,edem) = 1;
eq=eq+1;

%-----------------------
% [14] MARKUP shocks (Rescaled)
GAM0(eq,vmupR)      = 1;
GAM1(eq,vmupR)      = rhozetaMU;
PSI(eq,e_mupR)      = 1;
eq=eq+1;

%-----------------------------
% [15] Potential output (NO markup-shocks)
GAM0(eq,vynat) = 1/(1-Cap_phi_g1)+(varphi+alpha)/(1-alpha);
GAM0(eq,vg)    = -1/(1-Cap_phi_g1);
GAM0(eq,va)    = Cap_phi_g1/(1-Cap_phi_g1);

GAM1(eq,vynat) = Cap_phi_g1/(1-Cap_phi_g1);
GAM1(eq,vg)    = -Cap_phi_g1/(1-Cap_phi_g1);
eq=eq+1;

if nst>1
    if ZLB==0
        GAM0(eq:eq+nst-1,ves1:ves1+nst-1)    = eye(nst);
        GAM1(eq:eq+nst-1,ves1:ves1+nst-1)    = HB;
        PSI(eq:eq+nst-1,sh_es1:sh_es1+nst-1) = eye(nst);
        eq=eq+nst;
    else
        GAM0(eq,ves1)    = 1;
        GAM1(eq,ves1)    = 1;
        PSI(eq,sh_es1)   = 1;
        eq=eq+1;
    end
end

%-----------------------
% [16] EXP ERROR output
GAM0(eq,vy)    = 1;
GAM1(eq,vEy)   = 1;
PPI(eq,1) = 1;
eq=eq+1;

%-----------------------
% [17] EXP ERROR inflation
GAM0(eq,vpi)    = 1;
GAM1(eq,vEpi)   = 1;
PPI(eq,2) = 1;
eq=eq+1;


%-----------------------
% [18] EXP ERROR long term return
GAM0(eq,vRLT)  = 1;
GAM1(eq,vERLT) = 1;
PPI(eq,3) = 1;


end





