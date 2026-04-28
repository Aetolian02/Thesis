% =========================================================================
% CODE FOR "Methods for Measuring Expectations and Uncertainty in
% Markov-Switching Models"
% The Journal of Econometrics, 2016, 190(1), pp. 79–99" 
% Authors: Francesco Bianchi 
% =========================================================================
function [AUX,OMEGAT,CSIT,CSIT_1] = mss_lom_nonhomog6_note_plus_dsge_RQ6_rev(bm,RQ,HC,aux)
% Function to compute law of motion of first and second moments of a MS
% process that allows for a MS constant and heteroskedasticity 
% Model: y_t=bm*y_t-1+RQ*eps_t
% bm:  Autoregressive coefficients, EQUATIONS ALONG ROWS, constant must be last
%      coeffient
% RQ:  Product of square root of covariance matrix of the shocks and matrix
%      mapping the shocks into the endog variables
% HC:  Transition matrix
% aux: If 2 -> Compute law of motion of squared first moments (AUX{1-3})
%      If 1 -> Compute auxiliary output used to construct OMEGAT and CSIT
%      If 3 -> 1+2


% The code assumes a common chain for coeff and cov matrix


% The code allows for 4 levels:
% Level 1: Law of motion for first and second moments (OMEGAT and CSIT) and
%          auto-second moments (CSIT_1)
%          This part works only when a constant is NOT treated as a
%          variable always equal to one
%          The auxiliary output (AUX1) contains all the components used to
%          compute OMEGAT and CSIT: OMEGA, LOMpsi, M (->OMEGAT),
%              CSI, VV, LOMK, LOMPSI, DAC, MM (->CSIT). These last three
%              are not computed if a constant is not included
% Level 2: Law of motion for qq (OMkOM), computed if aux==2
% Level 3: Law of motion for qq1 (IkOM), computed if aux==2
% Level 4: Laws of motion for recursive sum of "auto" second and squared first moments
%          (LOMDDTQ and LOMDDTqq)





% Set some parameters
m = size(HC,2);   % Number of regimes
n = size(bm,1); % Number of variables
k = size(RQ,2); % Number of shocks
l = floor(size(bm,2)/n); % Number of lags


if size(bm,1)>size(bm,2);ddddddddddddd;end
% if size(bm,3)~=m;ddddddddddddd;end
if mod(size(bm,2),n)~=0; constant=1;bc=squeeze(bm(:,end,:));bm=bm(:,1:end-1,:); else constant=0;end
if l~=1;aaaaaaaaaaa;end % Second moments have not been asjusted fot l
if n==1;if l~=1; constant=1;bc=squeeze(bm(:,end,:));bm=bm(:,1:end-1,:); end; end

nout = nargout;

AUX=[];

if sum(sum(HC))~=size(HC,1)
    disp('---------------------------------------------')
    disp('Columns of transition matrix have to sum to 1')
    disp('---------------------------------------------')
end



%----------------------------------------------------%
%------------ LEVEL 1--------------------------------%
%----------------------------------------------------%
if nout>3 % You will need CSI_1 for level 3
    [OMEGA,CSI,CSI_1] = mss_lom_homog_plus5_dsge(bm,HC);
else
    if nout<=2
    [OMEGA] = mss_lom_homog_plus5_dsge(bm,HC);    
    else
    [OMEGA,CSI] = mss_lom_homog_plus5_dsge(bm,HC);
    end
end


if constant==1
    M=zeros(n*m*l,m);
    for isc=1:m
        M((isc-1)*n+1:isc*n,isc)=bc(:,isc);
    end
    OMEGAT=[[OMEGA,M*HC];[zeros(m,m*n*l),HC]];

else
    OMEGAT=OMEGA;
    M=[];
end



%----------------------------------------------------%


if nout>2 % Compute CSIT augmenting CSI
    VV = zeros(n*n*m,m);
    vIk=reshape(eye(k),k^2,1);
    for reg=1:m
        sigt = kron(RQ(:,:,reg),RQ(:,:,reg));                
        VV((reg-1)*n*n+1:reg*n*n,reg)=sigt*vIk;
    end


    if constant==1
        bdDAC = zeros(n*n*m,n*m);
        cc=zeros(m*n*n,m);
        for reg=1:m
            bmc=bc(:,reg);
            bml=bm(:,:,reg);
            bdDAC((reg-1)*n*n+1:reg*n*n,(reg-1)*n+1:reg*n)=kron(bml,bmc)+kron(bmc,bml);
            cc((reg-1)*n*n+1:reg*n*n,reg)=kron(bmc,bmc);
        end
        MATpi=(VV+cc)*HC;
        Lpsi=kron(HC,eye(n));
        DAC=bdDAC*Lpsi;
%                 
        CSITU=[CSI,DAC,MATpi];
        % 12/13/2015
        CSITD=[zeros(size(OMEGAT,1),size(CSI,2)),OMEGAT];
        CSIT=[CSITU;CSITD];

    else
        MATpi=VV*HC;
        CSITU=[CSI,MATpi];
        CSITD=[zeros(size(HC,1),size(CSI,2)),HC];
        CSIT=[CSITU;CSITD];
    end


end

if nout>3 % AUTOCOVARIANCE MATRIX

    if constant==0
        CSIT_1=CSI_1;
    else % AUGMENT CSI_1 when you have a constant
        IkC=zeros(n*n*m,n*m);
        CkI=zeros(n*n*m,n*m);
        for reg=1:m
            bmc=bc(:,reg);
            IkC((reg-1)*n*n+1:reg*n*n,(reg-1)*n+1:reg*n)=kron(eye(n),bmc);
            CkI((reg-1)*n*n+1:reg*n*n,(reg-1)*n+1:reg*n)=kron(bmc,eye(n));        
        end
        IkC_Lpsi=IkC*Lpsi;
        CkI_Lpsi=CkI*Lpsi;

        CSIT_1{1}=[[CSI_1{1},IkC_Lpsi];[zeros(n*m,n*n*m),Lpsi]];
        CSIT_1{2}=[[CSI_1{2},CkI_Lpsi];[zeros(n*m,n*n*m),Lpsi]];
    end

end
%----------------------------------------------------%
if aux(1)==1 || aux(1)==3
    % First 3 values of AUX are used for law of motion sq first mom
    AUX{4}=OMEGA;AUX{5}=M;AUX{6}=HC;
    if nout>2
        AUX{7}=[];
        AUX{8}=CSI;AUX{9}=VV;AUX{10}=MATpi;
        if constant==1
            AUX{11}=cc;AUX{12}=DAC;
        end
    end
end
%----------------------------------------------------%
%----- END OF LEVEL 1--------------------------------%
%----------------------------------------------------%


if nout>2 && aux>=2
    %----------------------------------------------------%
    %------------ LEVEL 2 -------------------------------%
    %----------------------------------------------------%
    % Law of motion of squared first moments
    OMkOM=kron(OMEGAT,OMEGAT);
    AUX{1}=OMkOM;
    %----------------------------------------------------%
    %----- END OF LEVEL 2 -------------------------------%
    %----------------------------------------------------%
end


if nout>3 && aux>=2
    %----------------------------------------------------%
    %------------ LEVEL 3 -------------------------------%
    %----------------------------------------------------%
    % Law of motion of auto squared first moments
    OMkI=kron(OMEGAT,eye((n+constant)*m));AUX{2}=OMkI;
    IkOM=kron(eye((n+constant)*m),OMEGAT);AUX{3}=IkOM;
    %----------------------------------------------------%
    %----- END OF LEVEL 3 -------------------------------%
    %----------------------------------------------------%

end


