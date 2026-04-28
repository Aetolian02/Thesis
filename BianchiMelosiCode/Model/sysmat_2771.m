%==================================================================
% Replication codes for
% The Dire Effects of the Lack of Monetary and Fiscal Coordination, 
% by Francesco Bianchi and Leonardo Melosi
% Journal of Monetary Economics, Volume 104, pp. 1-22. 
%==================================================================

% This code completes the construction of the system matrices for state-space representation

if msel(15)==1 && RC(1)==1 %Model Spcification in msel(15) and RC(1)=1: model is solved successfully
    nst=msel(10); % Number of regimes
    %% Variable indices
    ii=1;
    vy    = ii;ii=ii+1;%names(ii,:)=(' vy  ');  [1] Output
    vpi   = ii;ii=ii+1;%names(ii,:)=(' vpi ');  [2] Inflation
    vR    = ii;ii=ii+1;%names(ii,:)=(' vr  ');  [3] FFR
    vb    = ii;ii=ii+1;%names(ii,:)=(' vb  ');  [4] Debt
    
    vg    = ii;ii=ii+1;%names(ii,:)=(' vd  ');  [5] g (transormation of G/Y)
    vtau  = ii;ii=ii+1;%names(ii,:)=(' vtau');  [6] Taxes
    vz    = ii;ii=ii+1;%names(ii,:)=(' vz  ');  [7] TFP
    ve    = ii;ii=ii+1;%names(ii,:)=(' vtr ');  [8] Short term comp of transfers
    
    vmupR   = ii;ii=ii+1;%names(ii,:)=(' vtr '); [9] Mark-up shock (rescaled)
    vtp     = ii;ii=ii+1;%names(ii,:)=(' vtr ');  [10]Term premium (res in b eq)
    ve_LT   = ii;ii=ii+1;%names(ii,:)=(' vtr ');  [11]Long term comp of transfers
    vynat   = ii;ii=ii+1;%names(ii,:)=(' vy ');    [12] Natural or potential output
    
    vdem  = ii;ii=ii+1;%      [12] Demand shock
    vQm   = ii;ii=ii+1;%     [13] Price long term bond
    vRLT  = ii;ii=ii+1;%     [14] Return long term bond
    jj=1;
    vEy   = ii;ii=ii+1; %   [16] Expected output
    vEpi  = ii;ii=ii+1;%    [17] Expected inflation
    
    vERLT = ii;ii=ii+1;%    [19] Expected return maturity bond
    
    clear ii
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extract constants in case you have MS shocks
    % IMPORTANT: Dummies for regime changes have to be at the end
    [nrT,~,nrP]=size(TT);% If nrP=1=>Only MS shocks=>Model solved with gensys=>Extract constants
    
    if nst>1
        if nrP==1 % No changes in structural parameters
            % Extract constants
            CC=zeros(size(TT,1)-nst,nst);
            % Define regime switch shocks
            ms_shs3=zeros(nst,nst,nst);
            for ireg_st=1:nst
                estart=zeros(nst,1);estart(ireg_st)=1;
                for ireg_end=1:nst
                    eend=zeros(nst,1);eend(ireg_end)=1;
                    ms_shs3(:,ireg_end,ireg_st)=eend-PP*estart;
                end
            end
            estart=zeros(nst,1);estart(1)=1;
            for inst=1:nst
                CC(:,inst)=TT(1:end-nst,end-nst+1:end)*estart+RR(1:end-nst,end-nst+1:end)*ms_shs3(:,inst,1);
            end
            % Eliminate last nst variables (dummies for regime in place)
            TT=TT(1:end-nst,1:end-nst);
            RR=RR(1:end-nst,1:end-nst);
            nrT=nrT-nst;
            snames=snames(1,1:end-nst);
        else % Changes in structural parameters
            nst=size(TT,3);
            CC=reshape(TT(1:end-1,end,:),nrT-1,nst);
            % Eliminate last variable (dummy for regime in place)
            TT=TT(1:end-1,1:end-1,:);
            RR=RR(1:end-1,1:end-1,:);
            nrT=nrT-1;
            snames=snames(1,1:end-1);  
        end
    else
        CC=zeros(size(TT,1),1);
    end
    
    %% Augment matrices T and R
    ncR=size(RR,2);
    nrT2=nrT+1; % Number of additional series for lagged values
    
    TT2=zeros(nrT2,nrT2,nrP);
    RR2=zeros(nrT2,ncR,nrP);
    CC2=zeros(nrT2,nst);
    TT2(1:nrT,1:nrT,:)=TT;
    RR2(1:nrT,1:ncR,:)=RR;
    CC2(1:nrT,:)=CC;
    
    ii=1;
    vy1=nrT+ii;snames{nrT+ii}='y1';ii=ii+1;
    snames=snames(:,1:nrT2);
    
    for ii=1:nrP
        TT2(vy1,vy,ii)=1;
    end
    TT=TT2;
    RR=RR2;
    CC=CC2;
    
    %% observation equations indices
    eq_ygr = 1;
    eq_pi  = 2;
    eq_r   = 3;
    eq_b   = 4;
    eq_g   = 5;
    eq_tau = 6;
    eq_tr  = 7;
    eq_rFFR=8;
    eq_Pm=9; 
    ny = 9;
    
    % ================================================
    ZZ=zeros(ny,nrT);
    ZZ(eq_ygr,vy,1)    =  1;
    ZZ(eq_ygr,vy1,1)   = -1;
    ZZ(eq_ygr,vz,1)    =  1;
    ZZ(eq_pi,vpi,1)    =  1;
    ZZ(eq_r,vR,1)      =  1;
    ZZ(eq_b,vb,1)      =  1;
    ZZ(eq_g,vg,1)      =  1;  
    ZZ(eq_tau,vtau,1)  =  1;
    ZZ(eq_tr,ve,1)     =  1;
    ZZ(eq_rFFR,vR,1)   =  1; 
    ZZ(eq_rFFR,vEpi,1)   =  -1;
    ZZ(eq_Pm,vQm,1)   =   1;
    
    pi_steady     =  p_ss(2,1);
    gam_st        =  p_ss(3,1);
    b_steady      =  p_ss(4,1);
    g_steady      =  p_ss(5,1);
    tau_steady    =  p_ss(6,1);
    ave_mat     = p_ss(15,1);
    
    beta = p_ss(14,:);
    rr_steady = (1+gam_st) / beta-1;
    tr_steady = (1-1/beta)*b_steady+tau_steady-(1-1/g_steady);
    beta1 = 1/beta;
    rhowood=(1-1/ave_mat)*beta1;
    
    DD = zeros(ny,1);
    DD(eq_ygr,:)  = gam_st;
    DD(eq_pi,:)   = pi_steady;
    DD(eq_r,:)    = (pi_steady + rr_steady);
    DD(eq_b,1)    = b_steady;
    DD(eq_g,1)    = log(g_steady);
    DD(eq_tau,1)  = tau_steady;
    DD(eq_tr,1)   = tr_steady;
    DD(eq_rFFR)   =rr_steady;
    DD(eq_Pm)=log(1/(exp(pi_steady + rr_steady)-rhowood));
    
    if msel(5)==-1% Add the output gap as observable
        ny=ny+1;
        DD(ny,1)   = 0;
        ZZ(ny,vy,1)    =  1;%ZZ(ny,vynat,1)    =  -1;
    end
    nselhor=0;
    if nrP>1
        Tx2=sum(abs(TT)>10^-10,3); %
        [~,indt2]=find(sum(Tx2,2)'~=nst);
            [OM,CSI] = mss_lom_homog_plus5_dsge(TT(indt2,indt2,:),PP);
            leig=max(abs(eig(CSI)));
            if leig>=1;
                RC=0;
                error('No MSS in sysmat_2771');
            end        
    end
end

%% Construct the Variance and Covariance Matrix of Structural Shocks
is=1;
shindex=1:9;nex=9;
QQ = zeros(nex,nex,nsc);
for i=1:nsc
    QQ(:,:,i)=diag(p_er(shindex,i).^2);
end
VV=[];
%% Construct the Variance and Covariance Matrix of Structural Shocks
% OBSERVATION ERRORS
oeselect=1:ny+nselhor;
oerr=o_er(oeselect,:);
HH=diag(oerr.^2);
