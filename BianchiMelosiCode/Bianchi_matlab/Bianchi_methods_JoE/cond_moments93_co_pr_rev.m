% =========================================================================
% CODE FOR "Methods for Measuring Expectations and Uncertainty in
% Markov-Switching Models"
% The Journal of Econometrics, 2016, 190(1), pp. 79–99" 
% Authors: Francesco Bianchi 
% =========================================================================
function [ ME, SD, VAR, COVM, SECM, CORM ] = cond_moments93_co_pr_rev( OMEGA, CSIT, m, Thor, co, yall, prall )

% This function computes the law of motion for variances and standard
% deviations of a MS-DSGE (this code allows for a MS constant) 
% assuming that the economy starts from the vector y and prob pr
% If y and pr are not provided, the code assumes that the economy is in
% steady state and computed the moments for each of the regimes.
% OMEGA: Law of motion for first moments
% CSIT : Law of motion for second moments
% m : Number of regimes
% Thor : Time horizion
% co   : Dummy variable for constant, when 1 -> constant
% yall : Data (or starting values in case of impulse responses)
% prall: Starting probabilities (1/0 when interested in "conditional IR")

% ----------------------------------------------------------------------- %
% All results for first moments are in deviations from steady state, so the
% code can be used to compute impulse responses

n=(size(OMEGA,1)/m)-co; % Number of endogenous variables
k2=(size(CSIT,1)-co*(n*m+m))/m-n*n;k=k2.^0.5; % Number of exog shocks

nit=nargin;
if nit>5
    ss=size(yall,1);
    nrep=size(yall,3);
else
    ss=1;
    yall=zeros(1,n);
    nrep=1;
end
    
if nit>6
    all=0;
else % Compute law of motion conditional on each regime
    all=1;
    prall=eye(m);
end

nout = max(nargout,1);


%------------------------------------------------------------------------%
% disp('--------------------------------------------------')
% disp('------ Compute law of motion for moments  --------')
% disp('--------------------------------------------------')
% First moments
OMEGA_t(:,:,1)=eye((n+co)*m);
for tt=2:Thor+1
    OMEGA_t(:,:,tt)=OMEGA*OMEGA_t(:,:,tt-1);
end

% Second moments
if nout>1
    n2=n^2;

    n2k2m=n2*m+co*n*m+m;
    selQ=repmat(eye(n2),1,m);

    %---------------------------------------------------------------------
    CSIT_t=zeros(n2k2m,n2k2m,Thor+1);
    CSIT_t_temp_1=eye(n2k2m);
    CSIT_t(:,:,1)=CSIT_t_temp_1;
    for tt=2:Thor+1
        
        CSIT_t_temp=CSIT*CSIT_t_temp_1;
        CSIT_t(:,:,tt)=CSIT_t_temp;
        CSIT_t_temp_1=CSIT_t_temp;
%         if mod(tt,10)==0;disp(tt/Thor);end
    end

end


%------------------------------------------------------------------------%
% First moments
selq=repmat(eye(n),1,m);
me=zeros(n,ss,Thor+1,all*(m-1)+1);
ME=zeros(n,ss,Thor+1,all*(m-1)+1,nrep);

if nout>1
COVM=zeros(n*n,ss,Thor+1,all*(m-1)+1,nrep);
SECM=zeros(n*n,ss,Thor+1,all*(m-1)+1,nrep);
else
SECM=[];COVM=[];
end
% keyboard
for irep=1:nrep
    y=yall(:,:,irep);
    

    if all==1 % Compute conditional moments conditional on all regimes
              % This is what cond_moments3 does
        yrep=repmat(y,m,1);
        prmo=kron(eye(m),ones(ss,1));
%         keyboard
        if nout>1
            [ qtilda0, Qtilda0 ] = construct_q0_Q0_rev( yrep, prmo, co, k );
        else
            [ qtilda0 ] = construct_q0_Q0_rev( yrep, prmo, co, k );    
        end
        for th=1:Thor+1
            xx=reshape(selq*OMEGA_t(1:n*m,:,th)*qtilda0,n,ss,m);
            me(:,:,th,:)=xx;
        end
    else
        prmo=prall(:,:,irep);
        if nout>1
            [ qtilda0, Qtilda0 ] = construct_q0_Q0_rev( y, prmo, co, k );
        else
            [ qtilda0 ] = construct_q0_Q0_rev( y, prmo, co, k );    
        end
        for th=1:Thor+1
            me(:,:,th)=selq*OMEGA_t(1:n*m,:,th)*qtilda0;
        end
    end
    ME(:,:,:,:,irep)=me;


    %------------------------------------------------------------------------%
    % Second moments
    if nout>1
        QtildaT=zeros(n2k2m,ss,Thor+1,(m-1)*all+1);
        secmt=zeros(n2,ss,Thor+1,(m-1)*all+1);
        
        if all==1
            QtildaT_temp=zeros(n2k2m,ss,Thor+1,all*(m-1)+1);
            secmt_temp=zeros(n2,ss,Thor+1,(m-1)*all+1);
            for th=1:Thor+1
                xx=CSIT_t(:,:,th)*Qtilda0;
                QtildaT_temp(:,:,th,:)=reshape(xx,n2k2m,ss,m);
                
                xx2=reshape(selQ*xx(1:m*n2,:),n2,ss,m);
                secmt_temp(:,:,th,:)=xx2;
                
             
            end
            QtildaT=QtildaT_temp;
            secmt=secmt_temp;
            

        else
            for th=1:Thor+1
                QtildaT(:,:,th)=CSIT_t(:,:,th)*Qtilda0;
                secmt(:,:,th)=selQ*QtildaT(1:m*n2,:,th);
            end
        end




        
        if nit>5 % If starting values are actual data var=secm-firm^2

            SQFM=NaN(size(secmt));
            for ij=1:(m-1)*all+1
                for indt=1:ss
                    meij=squeeze(me(:,indt,:,ij));meijt=meij';
                    for tt=1:Thor+1
    %                     sqfm(:,indt,tt,ij)=diag(meij(:,tt)*meijt(tt,:));
    %                     keyboard
                        SQFM(:,indt,tt,ij)=reshape((meij(:,tt)*meijt(tt,:)),n2,1);
                    end
                end
            end
    
            covm=secmt-SQFM;
        else
            if co==0 % First moments are equal to zero
                covm=secmt;
            else
                SQFM=NaN(size(secmt));
                for ij=1:(m-1)*all+1
                    for indt=1:ss
                        meij=squeeze(me(:,indt,:,ij));meijt=meij';
                        for tt=1:Thor+1
                            SQFM(:,indt,tt,ij)=reshape((meij(:,tt)*meijt(tt,:)),n2,1);
                        end
                    end
                end
                covm=secmt-SQFM;
            end

        end
        %-------------------------
        secm=secmt;
        %-------------------------
        SECM(:,:,:,:,irep)=secm;
        COVM(:,:,:,:,irep)=covm;
    end
    
end

if nout>1
    extr=1:n+1:n2;
    COVM(:,:,1,:)=0; % To fix numerical approximation for horizon 0
    VAR=(COVM(extr,:,:,:));
    SD=VAR.^0.5;
    if nout>5 
        CORM=COVM;
        ind1=kron(1:n,ones(1,n));
        ind2=kron(ones(1,n),1:n);
        for ir=1:n2
            CORM(ir,:,:,:)=COVM(ir,:,:,:)./(SD(ind1(ir),:,:,:).*SD(ind2(ir),:,:,:));
        end
        CORM=squeeze(CORM);
    end
    COVM=squeeze(COVM);
    SECM=squeeze(SECM);
    SD  =squeeze(SD);
    VAR =squeeze(VAR);
end
ME=squeeze(ME);





end

