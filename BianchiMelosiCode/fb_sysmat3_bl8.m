

function [TT,QQ,RR,HH,DD,ZZ,VV,RC,CC,MAT]=fb_sysmat3_bl8(p_st,p_ss,p_er,o_er,PP,msel,MAT)
% Solve a Markov-switching DSGE model 
% by Francesco Bianchi, francesco.bianchi@duke.edu
% Code: 3/15/2021

nsb = size(p_st,2);
nsc = size(p_er,2);
CC=[];
snames=[];


if msel(15)==1
    %     
    if size(p_st,2)~=size(PP,2)
        disp('Number of regimes and size of p_st are not consistent')
        save zzz_bl8
        dddddddddddd
    end
    % 
    if msel(13)==1
        
    elseif sum(msel(13)==[2,25,4,5,6,65,7,8])

    elseif sum(msel(13)==[3,33,34,35,-35]) 
        
      
        znew=MAT{7};%
        if (msel(25)==-1 || msel(25)==-3) && size(znew,1)>2
            %             
            sizePP=size(PP,1);
            siznew=size(znew,1);
            if msel(6)==23 || msel(6)==24;nvar2=12;
            elseif msel(6)==26000;nvar2=12;
            elseif msel(6)==21000;nvar2=10;
            elseif msel(6)==25001;nvar2=10;
            elseif msel(6)==25002;nvar2=12;
            elseif msel(6)==10002;
                if msel(1)==1;nvar2=30;
                elseif msel(1)==2 || msel(1)==3 || msel(1)==4;
                    if p_st(8)==0;nvar2=24;else
                        nvar2=27;
                    end
                end
            else nvar2=2*size(MAT{1},1);
            end
            
            siznew10=siznew/nvar2;
            %             
            if siznew10>sizePP
                znew=znew(1:sizePP*nvar2,:);
            elseif siznew10<sizePP
                misreg=sizePP-siznew10;
                if misreg<=siznew10
                    znew=[znew;znew(end-nvar2*misreg+1:end,:)];
                else
                    
                    znew=[znew;repmat(znew(end-nvar2+1:end),sizePP-1,1)];
                    znew=znew(1:sizePP*nvar2);
                end
            end
        end
        if msel(13)==-35
            znew=[];
        end

        [TT,RR,RC,znew,snames] = fb_dsgesolve3_new7(p_st,p_ss,PP,msel,znew);%

    elseif msel(13)==-2 || msel(13)==-4 || msel(13)==-5 || msel(13)==-6 || msel(13)==-7 || msel(13)==-8
       
        
    elseif msel(13)==333
    
    end
    
    
    
elseif msel(15)==2 % Load solution
    

    [TT,~,RR,~,DD,ZZ,~,CC]=from_MAT_to_T(MAT);
    RC=1;
 
    
end





if RC(1)==1
    dxx=num2str(msel(6));
    sysmat_all=strcat('sysmat_',dxx);
    eval(sysmat_all)
    if msel(13)==34 && msel(15)==1% Check has to be done ex post when MS constant is modeled with dummy variables
        [~,CSI] = mss_lom_homog_plus5_dsge(TT,PP);
        try
            largestEig=eigs(CSI,1);
        catch
            largestEig=eig(CSI);

        end

        if largestEig<1
            RC = 1; % MSV is MSS
        else
            RC = 0; % MSV is NOT MSS
        end
    end
else
    QQ=[];HH=[];
    DD=[];ZZ=[];
    VV=[];CC=[];
end



if msel(15)==1

    
    if RC(1)==1
        MAT=from_T_to_MAT(TT,QQ,RR,HH,DD,ZZ,znew,CC,snames);
        RC=1;

    else
        MAT=[];
        RC=RC(1);
    end

elseif msel(15)==2
    [TT,~,RR,~,DD,ZZ,znew,CC,snames]=from_MAT_to_T(MAT);
    MAT=from_T_to_MAT(TT,QQ,RR,HH,DD,ZZ,znew,CC,snames);
    
    RC=1;
    
end





