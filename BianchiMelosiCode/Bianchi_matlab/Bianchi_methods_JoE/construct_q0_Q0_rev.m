% =========================================================================
% CODE FOR "Methods for Measuring Expectations and Uncertainty in
% Markov-Switching Models"
% The Journal of Econometrics, 2016, 190(1), pp. 79–99" 
% Authors: Francesco Bianchi 
% =========================================================================
function [ qtilda0, Qtilda0 ] = construct_q0_Q0_rev( y, prmo, co, k )
% This function construct starting values for first and second moments
% y = data
% prmo = regime probabilities
% co = dummy for constant (if 1, we have a MS constant)
% k = number of exogenous shocks;

n=size(y,2);n2=n^2;
m=size(prmo,2);
ss=size(y,1);

yt=y';
prt=prmo';
% Inital moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First moments
Un=ones(n,1);
fmy=repmat(yt,m,1);
prpr=kron(prt,Un);
q0=fmy.*prpr;
if co==1
    qtilda0=[q0;prt];
else
    qtilda0=q0;    
end

if nargout>1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Second moments
    smy=zeros(n2,ss);
    for ii=1:ss
        smy(:,ii)=reshape(yt(:,ii)*y(ii,:),n2,1);
    end
    Un2=ones(n2,1);
    smy=repmat(smy,m,1);
    prpr2=kron(prt,Un2);
    Q0=smy.*prpr2;


    if co==1
        Qtilda0=[Q0;qtilda0];
    else
        Qtilda0=[Q0;prt];
    end
end

end

