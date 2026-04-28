% =========================================================================
% CODE FOR "Methods for Measuring Expectations and Uncertainty in
% Markov-Switching Models"
% The Journal of Econometrics, 2016, 190(1), pp. 79–99" 
% Authors: Francesco Bianchi 
% =========================================================================
function [OMEGA,CSI,CSI_1] = mss_lom_homog_plus5_dsge(bm,P)
% Function to compute law of motion of first and second moments of a MS
% process that allows for more than one lag
% The code returns the law of motion for an homogenous system


% bm:       VAR coefficients (constant is assumed to be the last variable,
%           dependent variables along columns)
            % Modified for DSGE: Variables along lines
% P:        Transition matrix, columns sum to 1




ns=size(P,2);
n=size(bm,1);
l=floor(size(bm,2)/n);if l~=1;dddddddddddd;end
c=size(bm,2)-l*n;
if c~=0;dddddddddddddddd;end

nout = max(nargout,1);

if ns==1
    bm1=bm(:,:,1); % Exlude the constant (if there is one)
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=bm1;
    if nout>1
        bbm1=kron(bm1,bm1);
        CSI=bbm1;
    end
    
    if nout==3
        bbm1=kron(eye(n),bm1);
        CSI_1{1}=bbm1;

        bbm1=kron(bm1,eye(n));
        CSI_1{2}=bbm1;
    end
elseif ns==2
    bm1=bm(:,:,1); % Exlude the constant (if there is one)
    bm2=bm(:,:,2);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2)*(kron(P,eye(n*l)));
    
    if nout>1
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        CSI=blockdiag(bbm1,bbm2)*(kron(P,eye((n*l)^2)));
    end
%     keyboard
    if nout==3
        bbm1=kron(eye(n*l),bm1); % 5/29/2012: replace eye(n) with eye(n*l)
        bbm2=kron(eye(n*l),bm2); % 5/29/2012: replace eye(n) with eye(n*l)
        CSI_1{1}=blockdiag(bbm1,bbm2)*(kron(P,eye((n*l)^2)));

        bbm1=kron(bm1,eye(n*l)); % 5/29/2012: replace eye(n) with eye(n*l)
        bbm2=kron(bm2,eye(n*l)); % 5/29/2012: replace eye(n) with eye(n*l)
        CSI_1{2}=blockdiag(bbm1,bbm2)*(kron(P,eye((n*l)^2)));
    end
elseif ns==3
    bm1=bm(:,:,1);
    bm2=bm(:,:,2);
    bm3=bm(:,:,3);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3)*(kron(P,eye(n*l)));
    if nout>1
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        bbm3=kron(bm3,bm3);
        CSI=blockdiag(bbm1,bbm2,bbm3)*(kron(P,eye((n*l)^2)));
    end
    if nout==3
        bbm1=kron(eye(n*l),bm1);
        bbm2=kron(eye(n*l),bm2);
        bbm3=kron(eye(n*l),bm3);
        CSI_1{1}=blockdiag(bbm1,bbm2,bbm3)*(kron(P,eye((n*l)^2)));

        bbm1=kron(bm1,eye(n*l));
        bbm2=kron(bm2,eye(n*l));
        bbm3=kron(bm3,eye(n*l));
        CSI_1{2}=blockdiag(bbm1,bbm2,bbm3)*(kron(P,eye((n*l)^2)));
    end
elseif ns==4
    bm1=bm(:,:,1);
    bm2=bm(:,:,2);
    bm3=bm(:,:,3);
    bm4=bm(:,:,4);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4)*(kron(P,eye(n*l)));
    
    if nout>1
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        bbm3=kron(bm3,bm3);
        bbm4=kron(bm4,bm4);
        % keyboard
        CSI=blockdiag(bbm1,bbm2,bbm3,bbm4)*(kron(P,eye((n*l)^2)));
    end
    if nout==3
        bbm1=kron(eye(n*l),bm1);
        bbm2=kron(eye(n*l),bm2);
        bbm3=kron(eye(n*l),bm3);
        bbm4=kron(eye(n*l),bm4);
        CSI_1{1}=blockdiag(bbm1,bbm2,bbm3,bbm4)*(kron(P,eye((n*l)^2)));

        bbm1=kron(bm1,eye(n*l));
        bbm2=kron(bm2,eye(n*l));
        bbm3=kron(bm3,eye(n*l));
        bbm4=kron(bm4,eye(n*l));
        CSI_1{2}=blockdiag(bbm1,bbm2,bbm3,bbm4)*(kron(P,eye((n*l)^2)));
    end
elseif ns==5
    bm1=bm(:,:,1);
    bm2=bm(:,:,2);
    bm3=bm(:,:,3);
    bm4=bm(:,:,4);
    bm5=bm(:,:,5);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5)*(kron(P,eye(n*l)));
    
    if nout>1
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        bbm3=kron(bm3,bm3);
        bbm4=kron(bm4,bm4);
        bbm5=kron(bm5,bm5);
        % keyboard
        CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5)*(kron(P,eye((n*l)^2)));
    end
    if nout==3
        bbm1=kron(eye(n*l),bm1);
        bbm2=kron(eye(n*l),bm2);
        bbm3=kron(eye(n*l),bm3);
        bbm4=kron(eye(n*l),bm4);
        bbm5=kron(eye(n*l),bm5);
        CSI_1{1}=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5)*(kron(P,eye((n*l)^2)));

        bbm1=kron(bm1,eye(n*l));
        bbm2=kron(bm2,eye(n*l));
        bbm3=kron(bm3,eye(n*l));
        bbm4=kron(bm4,eye(n*l));
        bbm5=kron(bm5,eye(n*l));
        CSI_1{2}=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5)*(kron(P,eye((n*l)^2)));
    end
elseif ns==6
    bm1=bm(:,:,1);bm2=bm(:,:,2);
    bm3=bm(:,:,3);bm4=bm(:,:,4);
    bm5=bm(:,:,5);bm6=bm(:,:,6);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm6=[bm6;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5,bm6)*(kron(P,eye(n*l)));
    
    if nout>1      
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        bbm3=kron(bm3,bm3);
        bbm4=kron(bm4,bm4);
        bbm5=kron(bm5,bm5);
        bbm6=kron(bm6,bm6);
        % keyboard
        CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6)*(kron(P,eye((n*l)^2)));
    end
    if nout==3
        bbm1=kron(eye(n*l),bm1);
        bbm2=kron(eye(n*l),bm2);
        bbm3=kron(eye(n*l),bm3);
        bbm4=kron(eye(n*l),bm4);
        bbm5=kron(eye(n*l),bm5);
        bbm6=kron(eye(n*l),bm6);
        CSI_1{1}=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6)*(kron(P,eye((n*l)^2)));

        bbm1=kron(bm1,eye(n*l));
        bbm2=kron(bm2,eye(n*l));
        bbm3=kron(bm3,eye(n*l));
        bbm4=kron(bm4,eye(n*l));
        bbm5=kron(bm5,eye(n*l));
        bbm6=kron(bm6,eye(n*l));
        CSI_1{2}=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6)*(kron(P,eye((n*l)^2)));
    end
elseif ns==7
    bm1=bm(:,:,1);bm2=bm(:,:,2);
    bm3=bm(:,:,3);bm4=bm(:,:,4);
    bm5=bm(:,:,5);bm6=bm(:,:,6);bm7=bm(:,:,7);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm6=[bm6;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm7=[bm7;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end   
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5,bm6,bm7)*(kron(P,eye(n*l)));
    
    if nout>1
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        bbm3=kron(bm3,bm3);
        bbm4=kron(bm4,bm4);
        bbm5=kron(bm5,bm5);
        bbm6=kron(bm6,bm6);
        bbm7=kron(bm7,bm7);
        % keyboard
        CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6,bbm7)*(kron(P,eye((n*l)^2)));
    end
elseif ns==8
    bm1=bm(:,:,1);bm2=bm(:,:,2);
    bm3=bm(:,:,3);bm4=bm(:,:,4);
    bm5=bm(:,:,5);bm6=bm(:,:,6);
    bm7=bm(:,:,7);bm8=bm(:,:,8);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm6=[bm6;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm7=[bm7;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm8=[bm8;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5,bm6,bm7,bm8)*(kron(P,eye(n*l)));
    
    if nout>1
        bbm1=kron(bm1,bm1);
        bbm2=kron(bm2,bm2);
        bbm3=kron(bm3,bm3);
        bbm4=kron(bm4,bm4);
        bbm5=kron(bm5,bm5);
        bbm6=kron(bm6,bm6);
        bbm7=kron(bm7,bm7);
        bbm8=kron(bm8,bm8);
        % keyboard
        CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6,bbm7,bbm8)*(kron(P,eye((n*l)^2)));
    end
elseif ns==9
    bm1=bm(:,:,1);bm2=bm(:,:,2);
    bm3=bm(:,:,3);bm4=bm(:,:,4);
    bm5=bm(:,:,5);bm6=bm(:,:,6);
    bm7=bm(:,:,7);bm8=bm(:,:,8);bm9=bm(:,:,9);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm6=[bm6;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm7=[bm7;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm8=[bm8;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm9=[bm9;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5,bm6,bm7,bm8,bm9)*(kron(P,eye(n*l)));
    if nout>1
    bbm1=kron(bm1,bm1);
    bbm2=kron(bm2,bm2);
    bbm3=kron(bm3,bm3);
    bbm4=kron(bm4,bm4);
    bbm5=kron(bm5,bm5);
    bbm6=kron(bm6,bm6);
    bbm7=kron(bm7,bm7);
    bbm8=kron(bm8,bm8);
    bbm9=kron(bm9,bm9);
    % keyboard
    CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6,bbm7,bbm8,bbm9)*(kron(P,eye((n*l)^2)));
    end
elseif ns==10
    bm1=bm(:,:,1);bm2=bm(:,:,2);
    bm3=bm(:,:,3);bm4=bm(:,:,4);
    bm5=bm(:,:,5);bm6=bm(:,:,6);
    bm7=bm(:,:,7);bm8=bm(:,:,8);
    bm9=bm(:,:,9);bm10=bm(:,:,10);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm6=[bm6;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm7=[bm7;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm8=[bm8;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm9=[bm9;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm10=[bm10;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5,bm6,bm7,bm8,bm9,bm10)*(kron(P,eye(n*l)));
    if nout>1
    bbm1=kron(bm1,bm1);
    bbm2=kron(bm2,bm2);
    bbm3=kron(bm3,bm3);
    bbm4=kron(bm4,bm4);
    bbm5=kron(bm5,bm5);
    bbm6=kron(bm6,bm6);
    bbm7=kron(bm7,bm7);
    bbm8=kron(bm8,bm8);
    bbm9=kron(bm9,bm9);
    bbm10=kron(bm10,bm10);
    % keyboard
    CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6,bbm7,bbm8,bbm9,bbm10)*(kron(P,eye((n*l)^2)));
    end
elseif ns==11
    bm1=bm(:,:,1);bm2=bm(:,:,2);
    bm3=bm(:,:,3);bm4=bm(:,:,4);
    bm5=bm(:,:,5);bm6=bm(:,:,6);
    bm7=bm(:,:,7);bm8=bm(:,:,8);
    bm9=bm(:,:,9);bm10=bm(:,:,10);
    bm11=bm(:,:,11);
    if l>1
        bm1=[bm1;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm2=[bm2;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm3=[bm3;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm4=[bm4;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm5=[bm5;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm6=[bm6;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm7=[bm7;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm8=[bm8;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm9=[bm9;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm10=[bm10;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        bm11=[bm11;[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
    end
    OMEGA=blockdiag(bm1,bm2,bm3,bm4,bm5,bm6,bm7,bm8,bm9,bm10,bm11)*(kron(P,eye(n*l)));
    if nout>1
    bbm1=kron(bm1,bm1);bbm2=kron(bm2,bm2);
    bbm3=kron(bm3,bm3);bbm4=kron(bm4,bm4);
    bbm5=kron(bm5,bm5);bbm6=kron(bm6,bm6);
    bbm7=kron(bm7,bm7);bbm8=kron(bm8,bm8);
    bbm9=kron(bm9,bm9);bbm10=kron(bm10,bm10);
    bbm11=kron(bm11,bm11);
    % keyboard
    CSI=blockdiag(bbm1,bbm2,bbm3,bbm4,bbm5,bbm6,bbm7,bbm8,bbm9,bbm10,bbm11)*(kron(P,eye((n*l)^2)));
    end
elseif ns>=12
   
    OMEGAt=zeros(n*ns*l,n*ns*l);
    if nout>1
        CSIt=zeros(ns*(n*l)*(n*l),ns*(n*l)*(n*l));
    end
    
    for ins=1:ns
        if l==1
            bmt=bm(:,:,ins);
        else 
            check that this is correct using smaller number of regimes
            bmt=[bm(:,:,ins);[eye(n*(l-1),n*(l-1)),zeros(n*(l-1),n)]];
        end 
        
        st=(ins-1)*n*l+1;
        en=ins*n*l;
        OMEGAt(st:en,st:en)=bmt;
        if nout>1
            st=(ins-1)*(n*l)*(n*l)+1;
            en=ins*(n*l)*(n*l);
            bbmt=kron(bmt,bmt);
            CSIt(st:en,st:en)=bbmt;
        end
    end
    OMEGA=OMEGAt*(kron(P,eye(n*l)));
    if nout>1
    CSI=CSIt*(kron(P,eye((n*l)^2)));
    end
    
    
end












