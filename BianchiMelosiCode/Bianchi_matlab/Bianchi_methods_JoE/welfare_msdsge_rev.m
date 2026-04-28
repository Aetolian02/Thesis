% =========================================================================
% CODE FOR "Methods for Measuring Expectations and Uncertainty in
% Markov-Switching Models"
% The Journal of Econometrics, 2016, 190(1), pp. 79–99" 
% Authors: Francesco Bianchi 
% =========================================================================
function [ welfarelossT ] = welfare_msdsge_rev( OMEGA, CSIT, m, co, p_st, p_ss,oparams,posvar,yall, prall )

% This function computes welfare losses for a small size MS-DSGE 
% assuming that the economy starts from the vector y and prob pr
% If y and pr are not provided, the code assumes that the economy is in
% steady state and computes the welfare losses for each of the regimes.
% OMEGA: Law of motion for first moments
% CSIT : Law of motion for second moments
% m    : Number of regimes
% Thor : Time horizion
% co   : Dummy variable for constant, when 1 -> constant
% yall : Data (or starting values in case of impulse responses)
% prall: Starting probabilities (1/0 when interested in "conditional IR")

kappa=p_st(2,1);
beta=p_ss(14);
eps=oparams(7);
epi = posvar(1);
ey  = posvar(2);

n=(size(OMEGA,1)/m)-co; % Number of endogenous variables
k2=(size(CSIT,1)-co*(n*m+m))/m-n*n;k=k2.^0.5; % Number of exog shocks


evec = zeros(1,n);
evec(1,epi) = 1;
evec(1,ey)  = kappa/eps;

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


n2=n^2;
k2=size(CSIT,1)/m-n2-co*(n+1);%k=k2.^0.5;
selQ=repmat(eye(n2),1,m);


for irep=1:nrep
    y=yall(:,:,irep);
    

    if all==1 % Compute conditional moments conditional on all regimes
              % This is what cond_moments3 does
        yrep=repmat(y,m,1);
        prmo=kron(eye(m),ones(ss,1));

        [ ~, Qtilda0 ] = construct_q0_Q0_rev( yrep, prmo, co, k );

    else
        prmo=prall(:,:,irep);

        [ ~, Qtilda0 ] = construct_q0_Q0_rev( y, prmo, co, k );

    end



    %------------------------------------------------------------------------%
    % Second moments
    extr=1:n+1:n2;
    % This applies if only future squared deviations matter
    % COMP = (eye(size(CSIT,1))-beta*CSIT)\(beta*CSIT);
    % This applies if welfare depends on current deviation from SS
    % (CORRECT ONE): If you do not start from SS the two do not return the
    % same result
    COMP = (eye(size(CSIT,1))-beta*CSIT)\eye(size(CSIT,1));  
    Wtilda = [selQ,zeros(n2, m*k2+co*(m*n+m))];
    discsumQ = Wtilda*COMP*Qtilda0;
    welfarelossT = -evec*discsumQ(extr,:);
    
end

end

