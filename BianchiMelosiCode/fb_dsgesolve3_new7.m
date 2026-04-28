% ********************************************************************
% **   Solve MS-DSGE model using method by FWZ
% ********************************************************************

function [T1,T0,RC,znew,snames] = fb_dsgesolve3_new7(p_st,p_ss,PP,msel,varargin)
% This code solves a MS-DGSE using the NEW solution method of FWZ
% Code by Francesco Bianchi, francesco.bianchi@duke.edu

% Model in this form:
% A(s(t) x(t) = B(s(t) x(t-1) + Psi(s(t)) epsilon(t) + Pi eta(t)

% Solution method assumes invertible A matrix

% The solution is of the form
%
%   x(t) = V{s(t)}*F1{s(t)}*x(t-1) + V{s(t)}*G1{s(t)}*epsilon(t)
%
%  eta(t) = F2{s(t)}*x(t-1) + G2{s(t)}*epsilont(t)



% Dummy variable controlling the way solution method is initialized
sims=0;

nunu = num2str(msel(6)) ;
fcn=strcat('gammas_',nunu); 

h=size(p_st,2); % number of regimes
for i = 1:h

    [Afwz_cell{i},Bfwz_cell{i},gPsifwz_cell{i},PPI(:,:,i),snames] = feval(fcn,p_st(:,i),p_ss,msel);
end


l=size(PPI,2); % number of forward looking variables
maxit = 100; 
smallval = 1.0e-10;  %Convergence criterion: sqrt(machine epsilon).
ngensys = size(Afwz_cell{1},1);
nsh = size(gPsifwz_cell{1},2);





if min(size(varargin{1}))==0
    Is=eye(l);
    for ij=1:h


        [v111]=gensys_z_4(Afwz_cell{ij},Bfwz_cell{ij},PPI(:,:,ij),sims);
        v111=real(v111);
        [Qf{ij},Rf{ij}] = qr(Afwz_cell{ij}');
        Ui=[inv(Rf{ij}(1:ngensys-l,1:ngensys-l)'),zeros(ngensys-l,l)];
        U{ij}=Qf{ij}*[Ui;[zeros(l,ngensys-l),Is]];       
        UiViG=U{ij}\v111'; % General AiVi 
        UiVi=UiViG/UiViG(1:ngensys-l,:); % "Normalized" AiVi
        aa2=reshape(UiVi(end-l+1:end,:),l*(ngensys-l),1); % Reshape -Zi
        z00((ij-1)*l*(ngensys-l)+1:ij*l*(ngensys-l),1)=-aa2;
        
    end
else
    z00=varargin{1};
%     keyboard
end
% keyboard
if h>10 % Transition matrix likely to be sparse 
[F1fwz_cell, ~, G1fwz_cell, ~, Vfwz_cell, flag_err, znew] = fwz_msv_msre_singular_z3_sparse(PP', Afwz_cell, Bfwz_cell, gPsifwz_cell, l, z00, maxit,smallval); %flag_err
else
[F1fwz_cell, ~, G1fwz_cell, ~, Vfwz_cell, flag_err, znew] = fwz_msv_msre_singular_z3(PP', Afwz_cell, Bfwz_cell, gPsifwz_cell, l, z00, maxit,smallval); %flag_err    
end


if flag_err~=0 && min(size(varargin{1}))~=0
    Is=eye(l);
    for ij=1:h

        [v111]=gensys_z_4(Afwz_cell{ij},Bfwz_cell{ij},PPI(:,:,ij),sims);
        v111=real(v111);
        [Qf{ij},Rf{ij}] = qr(Afwz_cell{ij}');
        Ui=[inv(Rf{ij}(1:ngensys-l,1:ngensys-l)'),zeros(ngensys-l,l)];
        U{ij}=Qf{ij}*[Ui;[zeros(l,ngensys-l),Is]];       
        UiViG=U{ij}\v111'; % General AiVi 
        UiVi=UiViG/UiViG(1:ngensys-l,:); % "Normalized" AiVi
        aa2=reshape(UiVi(end-l+1:end,:),l*(ngensys-l),1); % Reshape -Zi
        z00((ij-1)*l*(ngensys-l)+1:ij*l*(ngensys-l),1)=-aa2;

    end

    if h>10 % Transition matrix likely to be sparse 
    [F1fwz_cell, ~, G1fwz_cell, ~, Vfwz_cell, flag_err, znew] = fwz_msv_msre_singular_z3_sparse(PP', Afwz_cell, Bfwz_cell, gPsifwz_cell, l, z00, maxit,smallval); %flag_err
    else
    [F1fwz_cell, ~, G1fwz_cell, ~, Vfwz_cell, flag_err, znew] = fwz_msv_msre_singular_z3(PP', Afwz_cell, Bfwz_cell, gPsifwz_cell, l, z00, maxit,smallval); %flag_err    
    end
end


if (flag_err)

    
    RC = -1; % NO MSV
    T1=0;
    T0=0;
%     
else
    

   T1=zeros(ngensys,ngensys,h);
   T0=zeros(ngensys,nsh,h);
   for si=1:h

      T1(:,:,si) = Vfwz_cell{si} * F1fwz_cell{si};
   end   
   for si=1:h

      T0(:,:,si)=Vfwz_cell{si} * G1fwz_cell{si};
   end
   
    if msel(13)==33
        %====== Checking stationarity of the solution. =======
        %   Checking the stationarity of regime-switching models. 

        nfwz = size(Afwz_cell{1},1);
        stackblkdiag = zeros(h*nfwz^2);
        for si=1:h
          stackblkdiag(((si-1)*nfwz^2+1):(si*nfwz^2),((si-1)*nfwz^2+1):(si*nfwz^2)) = kron(T1(:,:,si),T1(:,:,si));

        end 
        stackmat = kron(PP,eye(nfwz*nfwz)) * stackblkdiag; 
        


        if max(abs(eig(stackmat)))<1
            RC = 1; % MSV is MSS
        else 
            RC = 0; % MSV is NOT MSS
        end
    else
        RC=1;
    end
   
end

  





