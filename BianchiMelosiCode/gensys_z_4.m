function [z1,z2]=gensys_z_4(g0,g1,pi,sims,div)


%-------------------------------------------------------------
% Code by Francesco Bianchi based on gensys by Chris Sims to initialize FWZ solution method
% We impose the partition of the Z matrix in the QZ decomposition
% francesco.bianchi@duke.edu
%-------------------------------------------------------------
% function [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi,div)
% System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Christopher A. Sims
% Corrected 10/28/96 by CAS
% Modified by Francesco Bianchi to initialize FWZ
%---------------------------------------------------------------

if sims==1
%     eu=[0;0];
    realsmall=1e-7;
    fixdiv=(nargin==5);
    n=size(g0,1);
    [a b q z v]=qz(g0,g1);%keyboard
    if ~fixdiv, div=1.01; end
    nunstab=0;
    zxz=0;
    for i=1:n
    % ------------------div calc------------
       if ~fixdiv
          if abs(a(i,i)) > 0
             divhat=abs(b(i,i))/abs(a(i,i));
             if 1+realsmall<divhat & divhat<=div
                div=.5*(1+divhat);
             end
          end
       end
    % ----------------------------------------
       nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
       if abs(a(i,i))<realsmall & abs(b(i,i))<realsmall
          zxz=1;
       end
    end
    if ~zxz
       [a b q z]=qzdiv(div,a,b,q,z);
    end
else
    n=size(g0,1);
    [a b q z v]=qz(g0,g1); % QZ decomposition
    % ---------------------------------------------------------------
    % Reorder QZ decomposition, any of these two would work
    % div=1.01;[a b q z]=qzdiv(div,a,b,q,z);
    [a,b,q,z] = ordqz(a,b,q,z,'udo');
    
end
%----------------------------------------------------------------
nunstab=size(pi,2); % We impose number of nunstab roots
z1=z(:,1:n-nunstab)';
z2=z(:,n-nunstab+1:n)';
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%