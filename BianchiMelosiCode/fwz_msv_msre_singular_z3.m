function [F1, F2, G1, G2, V, err, z] = fwz_msv_msre_singular_z3(P, A, B, Psi, s, z, max_count, tol)
% Code by Francesco Bianchi to implement Farmer-Waggoner-Zha solution
% method for the case of a singular A matrix
% Code: 3/15/2021
%[ F1 F2 G1 G2 V err ] = fwz_msv_msre(P, A, B, Psi, s, x, max_count, tol)

% Computes MSV solution of 
%
%   A(s(t) x(t) = B(s(t) x(t-1) + Psi(s(t)) epsilon(t) + Pi eta(t)
%
% using Newton's Method.  Assumes that Pi' = [zeros(s,n-s)  eye(s)] and
% that A{i} is invertible.  P is the transition matrix and P(i,j) is the
% probability that s(t+1)=j given that s(t)=i.  Note that the rows of P
% must sum to one.  x is the initial value and if not passed is set to
% zero.  max_count is the maximum number of iterations on Newton's method
% before failing and tol is the convergence criterion.
%
% The solution is of the form
%
%   x(t) = V{s(t)}*F1{s(t)}*x(t-1) + V{s(t)}*G1{s(t)}*epsilon(t)
%
%  eta(t) = F2{s(t)}*x(t-1) + G2{s(t)}*epsilont(t)
%
% err is set to 0 upon sucess.  A positive value is the number of
% iterations before the method was terminated.
%


h=size(P,1);
n=size(A{1},1);
r=size(Psi{1},2);

if (nargin <= 7) || (tol <= 0)
    tol=1e-5;
end

if (nargin <= 6) || (max_count <= 0)
    max_count=1000;
end

if nargin <= 5
    x=zeros(h*s*(n-s),1);
end
f=zeros(h*s*(n-s),1);


Is=eye(s);
I=eye(n-s);

% The following code would be used if Pi were passed instead of the
% assumption that Pi' = [zeros(s,n-s)  eye(s)].
% for i=1:h
%     [Q,R] = qr(Pi{i});
%     W=[zeros(n-s,s)  I; inv(R(1:s,1:s)); zeros(s,n-s)]*Q';
%     A{i}=W*A{i};
%     B{i}=W*B{i};
%     Psi{i}=W*Psi{i};
% end

% save zzzz1
U=cell(h,1);
Qf=cell(h,1);
Rf=cell(h,1);
for i=1:h
%     U2{i}=inv(A{i});
    [Qf{i},Rf{i}] = qr(A{i}');
    Ui=[inv(Rf{i}(1:n-s,1:n-s)'),zeros(n-s,s)];
    U{i}=Qf{i}*[Ui;[zeros(s,n-s),Is]];
    AiUi = A{i}*U{i};
    C1{i} = AiUi(n-s+1:end,1:n-s);
    C2{i} = AiUi(n-s+1:end,n-s+1:end);
%     [Qf{i},Rf{i}] = qr(A{i});
    
end
% keyboard
C=cell(h,h);
for i=1:h
    for j=1:h
        C{i,j}=P(i,j)*B{j}*U{i};
    end
end

cont=true;
count=1;
D=zeros(h*s*(n-s),h*s*(n-s));
X=cell(h,1);
Z=cell(h,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you provide x
% z=[];
% for i=1:h
%     X{i}=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
%     Z{i}=pinv(C2{i})*(X{i}+C1{i});
% %     keyboard
%     z=[z;reshape(Z{i},size(Z{i},1)*size(Z{i},2),1)];
% %     keyboard
% %     disp('change this')
% %     Z{i}=zeros(size(Z{i}));
% %     X{i}=C2{i}*Z{i}-C1{i};
% %     disp('change this')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you provide z

for i=1:h
    Z{i}=reshape(z((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
    X{i}=C2{i}*Z{i}-C1{i};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard
% keyboard
while cont
%     for i=1:h
%         
%         for j=1:h
% %             W1=C{i,j}*[I; -X{i}];
%             W1=C{i,j}*[I; -Z{i}];
%             W2=W1(1:n-s,:);
% %             D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = kron(W2',Is);
%             D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = kron(W2',C2{j});
%             if i == j
%                 W1=zeros(s,n);
%                 for k=1:h
%                     W1=W1+[X{k}  Is]*C{i,k};
% %                     X2=C1{k}*Z{k}+C2{k};
% %                     W1=W1+[X2  Is]*C{i,k};
%                 end
%               W2=-W1(:,n-s+1:end);
%               D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) + kron(I,W2);
%             end
%         end
%     end
    
    for i=1:h
        mf=zeros(s,n-s);
        for j=1:h
%             mf=mf+[X{j} Is]*C{i,j}*[I; -X{i}];
            mf=mf+[X{j} Is]*C{i,j}*[I; -Z{i}];
        end
        f((i-1)*s*(n-s)+1:i*s*(n-s))=reshape(mf,s*(n-s),1);
        
    end
    
%     norm(f)
% %     real(f)
% %     keyboard
%     y=D\f;
% %     y=pinv(D)*f;
% %     x=x - y;
%     z=z - y;
%     keyboard
    if (count > max_count) || (norm(f) < tol)
        cont=false;
    else
        for i=1:h
        
            for j=1:h
    %             W1=C{i,j}*[I; -X{i}];
                W1=C{i,j}*[I; -Z{i}];
                W2=W1(1:n-s,:);
    %             D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = kron(W2',Is);
                D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = kron(W2',C2{j});
                if i == j
                    W1=zeros(s,n);
                    for k=1:h
                        W1=W1+[X{k}  Is]*C{i,k};
    %                     X2=C1{k}*Z{k}+C2{k};
    %                     W1=W1+[X2  Is]*C{i,k};
                    end
                  W2=-W1(:,n-s+1:end);
                  D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) + kron(I,W2);
                end
            end
        end
        % Modified by Francesco Bianchi on 11/22/2010
        % to avoid computing Newton's update when initial values already
        % satisfy norm(f)<tol
        % This is necesssary to avoid problems when using saved LOM matrices in
        % optimization with respect to parameters that do not affect LOM

        y=D\f;
%         keyboard
        z=z - y;
        count=count+1;
        for i=1:h
%             X{i}=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
% keyboard
            Z{i}=reshape(z((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
            X{i}=C2{i}*Z{i}-C1{i};
%             X{i}=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
%             Z{i}=inv(C2{i})*(X{i}+C1{i});
        end
    end
%     MF(cont,:)=mf;
end

if (norm(f) < tol)
    err=0;
else
    err=1;
end

F1=cell(h,1);
F2=cell(h,1);
G1=cell(h,1);
G2=cell(h,1);
V=cell(h,1);
pi=[zeros(n-s,s); Is]; 
% keyboard
for i=1:h
    Z=reshape(z((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
%     X=C2{i}*Z{i}-C1{i};
%     X=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
    V{i}=U{i}*[I; -Z];
    W=[A{i}*V{i} pi];
    F=W\B{i};
    F1{i}=F(1:n-s,:);
    F2{i}=F(n-s+1:end,:);
    
    G=W\Psi{i};
    G1{i}=G(1:n-s,:);
    G2{i}=G(n-s+1:end,:);
end






% count









  