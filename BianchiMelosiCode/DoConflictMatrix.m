%==================================================================
% Replication codes for
% The Dire Effects of the Lack of Monetary and Fiscal Coordination, 
% by Francesco Bianchi and Leonardo Melosi
% Journal of Monetary Economics, Volume 104, pp. 1-22. 
%==================================================================

function [HBl]=DoConflictMatrix(HB,oparams)

p_LD=oparams(1); %estimated probability that we stay out of the ZLB /high state of demand persists
pfightM=oparams(3); %peristence of the conflict regime 3 (fiscal backing removed temporarily)
pfightF=oparams(3); %peristence of the conflict regime 4 (fiscal backing removed long-lasting)
pfightFM=oparams(3); %peristence of the conflict regime 5 (mixed expectations about post-conflict policy mix backing removed forever)
% end
pF_M=oparams(4); %.7 %probability that after the conflict policymakers will coordinate on the ML  
% Transition Matrix
HBl=zeros(10,10);

%Coordinated out of the zlb  regimes (Regimes 1 and 2)
HBl(1:2,1:2)=HB(1:2,1:2);
HBl(6:9,1)=(1-p_LD)*.25; HBl(6:9,2)=(1-p_LD)*.25;
HBl(10,1)=0; HBl(10,2)=0;


%Fight Regimes (Regimes 3 and 4)
HBl(3,3)=p_LD*pfightM;  HBl(1,3)=p_LD*(1-pfightM); %exit to monetary led for sure
HBl(4,4)=p_LD*pfightF;  HBl(2,4)=p_LD*(1-pfightF); %exit to fiscally led for sure
HBl(5,5)=p_LD*pfightFM;  HBl(1,5)=p_LD*(1-pfightFM)*pF_M; HBl(2,5)=p_LD*(1-pfightFM)*(1-pF_M); %pF_M prob to exit to the monetary led regime 
HBl(6:10,3)=(1-p_LD)*.2; HBl(6:10,4)=(1-p_LD)*.2; HBl(6:10,5)=(1-p_LD)*.2;

%ZLB Regimes (Regime 5 and 7)
HBl(6,6)=HB(3,3); HBl(2,6)=(1-HB(3,3));  %Regime 6: coordinated exit ZLB-> FL
HBl(7,7)=HB(3,3); HBl(4,7)=(1-HB(3,3));  %Regime 7: non-coordinated exit (expected FL)
HBl(8,8)=HB(3,3); HBl(1,8)=(1-HB(3,3));   %Regime 8: coordinated exit ZLB->ML
HBl(9,9)=HB(3,3); HBl(3,9)=(1-HB(3,3));  %Regime 9: non-coordinated exit (expected ML)
HBl(10,10)=HB(3,3); HBl(5,10)=(1-HB(3,3));  %Regime 10: non-coordinated exit mixed expectations

