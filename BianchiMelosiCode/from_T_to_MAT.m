function [ MAT ] = from_T_to_MAT(Tm,Qm,Rm,Hm,Dm,Zm,znew,varargin)
% by Francesco Bianchi
% Code: 3/15/2021

MAT={Tm,Qm,Rm,Hm,Dm,Zm,znew};
% 
if nargin>7
    MAT{8}=varargin{1};
    if nargin>8
        MAT{9}=varargin{2};
        if nargin>9
            MAT{10}=varargin{3};
            if nargin>10
                MAT{11}=varargin{4};
            end
        end
    end
end


end

