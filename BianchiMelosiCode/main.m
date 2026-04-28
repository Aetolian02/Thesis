%==================================================================
% Replication codes for
% The Dire Effects of the Lack of Monetary and Fiscal Coordination, 
% by Francesco Bianchi and Leonardo Melosi
% Journal of Monetary Economics, Volume 104, pp. 1-22. 
% FB: francesco.bianchi@duke.edu
% LM: leonardo.melosi@chi.frb.org
%==================================================================
clear
close all
clc

addpath Documents/Licencjat/Bianchi_code/Bianchi_matlab
addpath Documents/Licencjat/Bianchi_code/Bianchi_matlab/Bianchi_methods_JoE
addpath Documents/Licencjat/Bianchi_code/Model
addpath Documents/Licencjat/Bianchi_code/CondMoms

%==========================================================================
%% Editable Part Begin
%% Customizing IRF
b0=.50; %initial stock of debt (above steady state)
nirf=31; % length of the IRF simulation
xx=10; % the length of ZLB period in the simulations
ft=10; %Duration of conflict period
%% Graphic settings
set(0,'defaultlinelinewidth',1.3);
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',12);
set(0,'defaultlinelinewidth',2.2)

rrs = get(0,'screensize');
rrf = get(0,'defaultfigurepos');
m = min((rrs(3)-200)/rrf(3)/2, (rrs(4)-200)/rrf(4)/2);
defpos = [100, 100, 1300, 550];
set(0,'defaultfigurepos',defpos);
bw=0; %bw=1 black and white graphs for paper bw=0 for slides
dufig    = 1; % Set equal to 1 to save figures in the subfolder Figures

cond_moments=1; %if set to 1, it computes conditional moments and save them (time consuming); else upload the moments from pre-existing files

%% Editable Part End
%==========================================================================

%% Upload Parameters
p_st=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','p_st');
p_ss=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','p_ss');
p_er=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','p_er');
o_er=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','o_er');
oparams=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','OtherParams');
HB=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','TransitionMatrix');
msel=xlsread('Documents/Licencjat/Bianchi_code/Params/Params.xlsx','msel'); % This vector contains model options 
% msel(6):  Number to select the model (2771)
% msel(13): Solution method: 35-> (Farmer Waggoner Zha) or gensys in this
%               paper and check MSS
% msel(15): 1->Solve model
%           2->Only build cov matrix and mapping from model to observables

MATlc{7}=[];  

[HBlc]=DoConflictMatrix(HB,oparams);

%% Model Solution
[TT,QQ,RR,HH,DD,ZZ,VV,RC,CC,MATmo]=...
    fb_sysmat3_bl8(p_st,p_ss,p_er,o_er,HBlc,msel,MATlc);


%% Simulations Conditional On Regime Paths
ssirfAF=zeros(size(TT,1),nirf); yyirfAF=zeros(size(ZZ,1),nirf);
ssirfAF(4,1)=b0*4; yyirfAF=DD+ZZ*ssirfAF(:,1);
ssirfLCF=zeros(size(TT,1),nirf); yyirfLCF=zeros(size(ZZ,1),nirf);
ssirfLCF(4,1)=b0*4; yyirfLCF=DD+ZZ*ssirfLCF(:,1);
ssirfAM=zeros(size(TT,1),nirf); yyirfAM=zeros(size(ZZ,1),nirf);
ssirfAM(4,1)=b0*4; yyirfAM=DD+ZZ*ssirfAM(:,1);
ssirfLCM=zeros(size(TT,1),nirf); yyirfLCM=zeros(size(ZZ,1),nirf);
ssirfLCM(4,1)=b0*4; yyirfLCM=DD+ZZ*ssirfLCM(:,1);
ssirfLCFMM=zeros(size(TT,1),nirf); yyirfLCFMM=zeros(size(ZZ,1),nirf);
ssirfLCFMM(4,1)=b0*4; yyirfLCFMM=DD+ZZ*ssirfLCFMM(:,1);
ssirfLCFMF=zeros(size(TT,1),nirf); yyirfLCFMF=zeros(size(ZZ,1),nirf);
ssirfLCFMF(4,1)=b0*4; yyirfLCFMF=DD+ZZ*ssirfLCFMF(:,1);
for t=2:xx+1 %ZLB period
    ssirfAF(:,t)=CC(:,6)+TT(:,:,6)*ssirfAF(:,t-1);
    yyirfAF(:,t)=DD+ZZ*ssirfAF(:,t);
    ssirfLCF(:,t)=CC(:,7)+TT(:,:,7)*ssirfLCF(:,t-1);
    yyirfLCF(:,t)=DD+ZZ*ssirfLCF(:,t);
    ssirfAM(:,t)=CC(:,8)+TT(:,:,8)*ssirfAM(:,t-1);
    yyirfAM(:,t)=DD+ZZ*ssirfAM(:,t);
    ssirfLCM(:,t)=CC(:,9)+TT(:,:,9)*ssirfLCM(:,t-1);
    yyirfLCM(:,t)=DD+ZZ*ssirfLCM(:,t);
    ssirfLCFMM(:,t)=CC(:,10)+TT(:,:,10)*ssirfLCFMM(:,t-1);
    yyirfLCFMM(:,t)=DD+ZZ*ssirfLCFMM(:,t);
    ssirfLCFMF(:,t)=CC(:,10)+TT(:,:,10)*ssirfLCFMF(:,t-1);
    yyirfLCFMF(:,t)=DD+ZZ*ssirfLCFMF(:,t);
end
for t=xx+2:ft+xx+1 %Fight period
    ssirfAF(:,t)=CC(:,2)+TT(:,:,2)*ssirfAF(:,t-1);
    yyirfAF(:,t)=DD+ZZ*ssirfAF(:,t);
    ssirfLCF(:,t)=CC(:,4)+TT(:,:,4)*ssirfLCF(:,t-1); %fight
    yyirfLCF(:,t)=DD+ZZ*ssirfLCF(:,t);
    ssirfAM(:,t)=CC(:,1)+TT(:,:,1)*ssirfAM(:,t-1);
    yyirfAM(:,t)=DD+ZZ*ssirfAM(:,t);
    ssirfLCM(:,t)=CC(:,3)+TT(:,:,3)*ssirfLCM(:,t-1); %fight
    yyirfLCM(:,t)=DD+ZZ*ssirfLCM(:,t);
    ssirfLCFMM(:,t)=CC(:,5)+TT(:,:,5)*ssirfLCFMM(:,t-1);
    yyirfLCFMM(:,t)=DD+ZZ*ssirfLCFMM(:,t);
    ssirfLCFMF(:,t)=CC(:,5)+TT(:,:,5)*ssirfLCFMF(:,t-1);
    yyirfLCFMF(:,t)=DD+ZZ*ssirfLCFMF(:,t);
end
for t=ft+xx+2:nirf %Coordination
    ssirfAF(:,t)=CC(:,2)+TT(:,:,2)*ssirfAF(:,t-1);
    yyirfAF(:,t)=DD+ZZ*ssirfAF(:,t);
    ssirfLCF(:,t)=CC(:,2)+TT(:,:,2)*ssirfLCF(:,t-1); %fight
    yyirfLCF(:,t)=DD+ZZ*ssirfLCF(:,t);
    ssirfAM(:,t)=CC(:,1)+TT(:,:,1)*ssirfAM(:,t-1);
    yyirfAM(:,t)=DD+ZZ*ssirfAM(:,t);
    ssirfLCM(:,t)=CC(:,1)+TT(:,:,1)*ssirfLCM(:,t-1); %fight
    yyirfLCM(:,t)=DD+ZZ*ssirfLCM(:,t);
    ssirfLCFMM(:,t)=CC(:,1)+TT(:,:,1)*ssirfLCFMM(:,t-1);
    yyirfLCFMM(:,t)=DD+ZZ*ssirfLCFMM(:,t);
    ssirfLCFMF(:,t)=CC(:,2)+TT(:,:,2)*ssirfLCFMF(:,t-1);
    yyirfLCFMF(:,t)=DD+ZZ*ssirfLCFMF(:,t);
end

nsb=size(TT,3);
for ib=1:nsb
    TTH(:,:,ib)=[TT(:,:,ib),CC(:,ib)];
    RQH(:,:,ib)=RR(:,:,ib)*QQ.^0.5;
end

[~,OMEGAH,CSIH] = mss_lom_nonhomog6_note_plus_dsge_RQ6_rev(TTH,RQH,HBlc,0);


%% Computing Conditional Moments
%   LCF path
smbcLCF=zeros(size(ssirfLCF,2)-1,nsb);
smbcLCF(1:10,7)=ones(10,1); %LD (FA is expected to win after conflict)
smbcLCF(11:20,4)=ones(10,1); %conflict (LB is expected to win)
smbcLCF(21:30,2)=ones(10,1);

%LCM Path
smbcLCM=zeros(size(ssirfLCM,2)-1,nsb);
smbcLCM(1:10,9)=ones(10,1); %LD (FA is expected to win after conflict)
smbcLCM(11:20,3)=ones(10,1); %conflict (LB is expected to win)
smbcLCM(21:30,1)=ones(10,1);
%   AF path
smbcAF=zeros(size(ssirfLCF,2)-1,nsb);
smbcAF(1:10,6)=ones(10,1); %ZLB (LB is expecting no conflict)
smbcAF(11:30,2)=ones(20,1); %FL

%AM Path
smbcAM=zeros(size(ssirfLCF,2)-1,nsb);
smbcAM (1:10,8)=ones(10,1); %ZLB (FA is expecting no conflict)
smbcAM (11:30,1)=ones(20,1); %FL
Thor=40;

if cond_moments==1
    [ ~, SD_LCF, ~, ~, ~, ~ ] = cond_moments93_co_pr_rev( OMEGAH, CSIH, nsb*1, Thor, 1, ssirfLCF(:,2:end)', smbcLCF );
    [ ~, SD_LCM, ~, ~, ~, ~ ] = cond_moments93_co_pr_rev( OMEGAH, CSIH, nsb*1, Thor, 1, ssirfLCM(:,2:end)', smbcLCM );
    [ ~, SD_AF, ~, ~, ~, ~ ] = cond_moments93_co_pr_rev( OMEGAH, CSIH, nsb*1, Thor, 1, ssirfAF(:,2:end)', smbcAF );
    [ ~, SD_AM, ~, ~, ~, ~ ] = cond_moments93_co_pr_rev( OMEGAH, CSIH, nsb*1, Thor, 1, ssirfAM(:,2:end)', smbcAM );
    save CondMoms\SD SD_LCF SD_LCM SD_AF SD_AM
else
    load SD SD_LCF SD_LCM SD_AF SD_AM
end

%Computing Welfare under Baseline Parameterization
posvar=[2;12]; % set the position of the variables in the quadratic function of welfare
[ welfarelossLCF ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1,p_st, p_ss,oparams,posvar, ssirfLCF(:,2:end)', smbcLCF );
[ welfarelossLCM ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_ss,oparams,posvar, ssirfLCM(:,2:end)', smbcLCM );
[ welfarelossAF ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_ss,oparams,posvar, ssirfAF(:,2:end)', smbcAF );
[ welfarelossAM ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_ss,oparams,posvar, ssirfAM(:,2:end)', smbcAM );

% Computing welfare assuming lower beta 
p_sslb=p_ss; p_sslb(14)=oparams(6);
[ welfarelossLCFbeta ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_sslb,oparams,posvar, ssirfLCF(:,2:end)', smbcLCF );
[ welfarelossLCMbeta ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_sslb,oparams,posvar, ssirfLCM(:,2:end)', smbcLCM );
[ welfarelossAFbeta ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_sslb,oparams,posvar, ssirfAF(:,2:end)', smbcAF );
[ welfarelossAMbeta ] = welfare_msdsge_rev( OMEGAH, CSIH, nsb*1, 1, p_st, p_sslb,oparams,posvar, ssirfAM(:,2:end)', smbcAM );

yyirfAM(4,:)=yyirfAM(4,:)/4;
yyirfLCM(4,:)=yyirfLCM(4,:)/4;
yyirfAF(4,:)=yyirfAF(4,:)/4;
yyirfLCF(4,:)=yyirfLCF(4,:)/4;
yyirfLCFMM(4,:)=yyirfLCFMM(4,:)/4;
yyirfLCFMF(4,:)=yyirfLCFMF(4,:)/4;

%% Producing Figures 3 and 4

%% Figure 3 of the paper
figure('name','Non-Coordinated exit (->FL) - Observables','numbertitle','off')
subplot(2,3,1)
plot(0:nirf-1,ssirfAF(1,:)*100,'--k');
hold on
if bw==1
    plot(0:nirf-1,ssirfLCF(1,:)*100,'k');
    hold on
    patch([1 xx xx 1],[-5 -5 2.5 2.5 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-5 -5 2.5 2.5 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,ssirfLCF(1,:)*100);
    hold on
    patch([1 xx xx 1],[-5 -5 2.5 2.5 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-5 -5 2.5 2.5 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ssirfAF(1,:)*0,'-.r','LineWidth',1);
hold off
axis tight
title('Output Gap') 

subplot(2,3,2)
plot(0:nirf-1,yyirfAF(2,:)*400,'--k');
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(2,:)*400,'k');
    hold on
    patch([1 xx xx 1],[0 0 5 5],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[0 -0 5 5],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
else
    plot(0:nirf-1,yyirfLCF(2,:)*400);
    hold on
    patch([1 xx xx 1],[0 0 5 5],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[0 0 5 5],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
end
hold on
plot(0:nirf-1,ones(nirf)*DD(2)*400,'-.r','LineWidth',1);
hold off

axis tight
title('Inflation') 


subplot(2,3,3)
plot(0:nirf-1,yyirfAF(3,:)*400,'--k');
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(3,:)*400,'k');
    hold on
    patch([1 xx xx 1],[0 0 7.5 7.5],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[0 0 7.5 7.5 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,yyirfLCF(3,:)*400);
    hold on
    patch([1 xx xx 1],[0 0 7.5 7.5],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[0 0 7.5 7.5 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(3)*400,'-.r','LineWidth',1);
hold off
axis([0 nirf-1 0 7.5])
title('FFR ') 


subplot(2,3,4)
plot(0:nirf-1,yyirfAF(4,:)*100,'--k');
hold on

if bw==1
    plot(0:nirf-1,yyirfLCF(4,:)*100,'k');
    hold on
    patch([1 xx xx 1],[75 75 120 120 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[75 75 120 120 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,yyirfLCF(4,:)*100);
    hold on
    patch([1 xx xx 1],[75 75 120 120 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[75 75 120 120 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(4)*25,'-.r','LineWidth',1);
hold off

axis([0 nirf-1 75 120])

title('Debt-to-GDP Ratio') 
h1=legend('LD->FL','LD->Conflict(FL)->FL');
set(h1,'Location','NorthWest','FontSize',10);

subplot(2,3,6)
plot(0:nirf-1,(yyirfAF(8,:))*400,'--k');
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(8,:)*400,'k');
    hold on
    patch([1 xx xx 1],[-1 -1 4 4],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-1 -1 4 4 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,(yyirfLCF(8,:))*400);
    hold on
    patch([1 xx xx 1],[-1 -1 4 4],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-1 -1 4 4],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(8)*400,'-.r','LineWidth',1);
hold off
axis([0 nirf-1 -1 4])
title('Real FFR ') 

subplot(2,3,5)
plot(0:nirf-1,(yyirfAF(9,:)),'--k');
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(9,:),'k');
    hold on
    patch([1 xx xx 1],[2 2 3 3],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[2 2 3 3 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,(yyirfLCF(9,:)));
    hold on
    patch([1 xx xx 1],[2.5 2.5 3 3],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[2.5 2.5 3 3 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(9),'-.r','LineWidth',1);
hold off
axis([0 nirf-1 2.5 3])
title('Price of Bonds ') 
if dufig==1
    savefigure_pdf('Figures\Figure3') ;
end

%% Figure 4

figure('name','Non-Coordinated exit (->ML)- Observables_2 ','numbertitle','off')
subplot(2,3,1)
plot(0:nirf-1,ssirfLCM(1,:)*100,'-ok','MarkerSize',3);
hold on
if bw==1
    plot(0:nirf-1,ssirfLCF(1,:)*100,'k');
    hold on
    patch([1 xx xx 1],[-5 -5 2.5 2.5 ],'k','EdgeColor','k','FaceAlpha',.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-5 -5 2.5 2.5 ],'k','EdgeColor','k','FaceAlpha',.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,ssirfLCF(1,:)*100);
    hold on
    patch([1 xx xx 1],[-5 -5 2.5 2.5 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-5 -5 2.5 2.5 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ssirfAF(1,:)*0,'-.r','LineWidth',1);
hold off
axis tight
title('Output Gap') 


subplot(2,3,2)
plot(0:nirf-1,yyirfLCM(2,:)*400,'o-k','MarkerSize',3);
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(2,:)*400,'k');
    hold on
    patch([1 xx xx 1],[-.5 -.5 5 5 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-0.5 -0.5 5 5 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
else
    plot(0:nirf-1,yyirfLCF(2,:)*400);
    hold on
    patch([1 xx xx 1],[-.5 -.5 5 5 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-0.5 -0.5 5 5 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
end
hold on
plot(0:nirf-1,ones(nirf)*DD(2)*400,'-.r','LineWidth',1);
hold off
axis([0 nirf-1 -.5 5])
title('Inflation') % (annualized in %)


subplot(2,3,3)
plot(0:nirf-1,yyirfLCM(3,:)*400,'o-k','MarkerSize',3);
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(3,:)*400,'k');
    hold on
    patch([1 xx xx 1],[0 0 7.5 7.5],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[0 0 7.5 7.5 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,yyirfLCF(3,:)*400);
    hold on
    patch([1 xx xx 1],[0 0 7.5 7.5],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[0 0 7.5 7.5 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(3)*400,'-.r','LineWidth',1);
hold off
axis([0 nirf-1 0 7.5])
title('FFR ') %(annualized in %)


subplot(2,3,4)
plot(0:nirf-1,yyirfLCM(4,:)*100,'o-k','MarkerSize',3);

hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(4,:)*100,'k');
    hold on
    patch([1 xx xx 1],[75 75 120 120 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[75 75 120 120 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,yyirfLCF(4,:)*100);
    hold on
    patch([1 xx xx 1],[75 75 120 120 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[75 75 120 120 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(4)*25,'-.r','LineWidth',1);
hold off
axis([0 nirf-1 75 120])

title('Debt-to-GDP Ratio') 
h1=legend('LD->Conflict(ML)->ML','LD->Conflict(FL) ->FL');
set(h1,'Location','SouthEast','FontSize',10);

subplot(2,3,6)
plot(0:nirf-1,yyirfLCM(8,:)*400,'o-k','MarkerSize',3);
hold on
if bw==1
    plot(0:nirf-1,yyirfLCF(8,:)*400,'k');
    hold on
    patch([1 xx xx 1],[-1 -1 4 4],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-1 -1 4 4 ],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,yyirfLCF(8,:)*400);
    hold on
    patch([1 xx xx 1],[-1 -1 4 4],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[-1 -1 4 4 ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(8)*400,'-.r','LineWidth',1);
hold off
axis([0 nirf-1 -1 4])
title('Real FFR ') 

subplot(2,3,5)
plot(0:nirf-1,yyirfLCM(9,:),'o-k','MarkerSize',3);
hold on


if bw==1
    plot(0:nirf-1,yyirfLCF(9,:),'k');
    hold on
    patch([1 xx xx 1],[2 2 3 3],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[2 2 3 3],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
else
    plot(0:nirf-1,(yyirfLCF(9,:)));
    hold on
    patch([1 xx xx 1],[2.5 2.5 3 3 ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[2.5 2.5 3 3  ],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
    hold on
end
plot(0:nirf-1,ones(nirf)*DD(9),'-.r','LineWidth',1);
hold off
axis([0 nirf-1 2.5 3])
title('Price of Bonds ') 
if dufig==1
    savefigure_pdf('Figures\Figure4') ;
end

%% The Model with the Shock-Specific Rule (Emergency Budget)
p_st_scsp=repmat(p_st(1:40,1),1,2);
p_st_scsp(3,2)=p_st(3,2)*0; p_st_scsp(8,2)=p_st(8,2)*0;
p_st_scsp(35,2)=p_st(35,6); %discrete demand shock

p_st_scsp(17,:)=0; 

HBss=zeros(2,2); HBss(1,1)=oparams(1); HBss(2,1)=1-oparams(1); HBss(1,2)=1-oparams(2); HBss(2,2)=oparams(2);
[GAM0scsp,GAM1scsp,PSIscsp,PPIscsp,~] = gammas_2771_shock_specific_Rules_4(p_st_scsp,p_ss,msel,HBss);
[TTshsp0,~,RRshsp0,~,~,~,~,eu,~]=gensys(GAM0scsp,GAM1scsp,zeros(size(GAM0scsp,1),1),PSIscsp,PPIscsp);
if sum(eu)~=2
    error('indeterminacy when solving the model with shock-specific policy rules')
end
%Remove the dummy variables
CCshsp(:,1)=RRshsp0(:,10); CCshsp(:,2)=RRshsp0(:,11);
%Resizing
CCshsp=CCshsp([1:23,26:size(CCshsp,1)],:);
TTshsp=TTshsp0([1:23,26:size(TTshsp0,1)],[1:23,26:size(TTshsp0,1)]);
RRshsp=RRshsp0([1:23,26:size(TTshsp0,1)],1:9);

ZZshsp=zeros(size(TTshsp)); ZZshsp(1,1:size(ZZ,2))=ZZ(1,:); ZZshsp(2:16,2:16)=eye(15);
ZZshsp(8,:)=ZZshsp(8,:)*0; ZZshsp(8,3)=1; ZZshsp(8,27)=-1;  %the real rate in the actual economy
ZZshsp(9,:)=ZZshsp(9,:)*0; ZZshsp(9,7)=1;   %the price of long-term bonds

DDshsp=[DD;zeros(size(TTshsp,1)-length(DD),1)]; DDshsp(10:13,1)=DD(1:4,1); DDshsp(12)=DDshsp(13);
QQshsp=QQ;
% Solve the  model where AM is announced (no fight)
p_st_lc4=p_st;
%Active regime has to be the same as in the shock-specific rule
p_st_lc4(1:40,1)=p_st_scsp(1:40,1); 

HBlc4=zeros(2,2); HBlc4(1,1)=oparams(1); HBlc4(2,1)=1-oparams(1);
HBlc4(2,2)=oparams(2); HBlc4(1,2)=1-oparams(2);

p_st_lc4=p_st_lc4(:,[1,7]);
%No FL even during the low-demand spell
p_st_lc4([3:5,8,10],2)=p_st_lc4([3:5,8,10],1);
p_st_lc4(1,2)=0;
[TT4,QQ4,RR4,HH4,DD4,ZZ4,VV4,RC4,CC4,MATmo4]=...
    fb_sysmat3_bl8(p_st_lc4,p_ss,p_er,o_er,HBlc4,msel,MATlc);


ssirfAM4=zeros(size(TT,1),nirf); yyirfAM4=zeros(size(ZZ,1),nirf);
ssirfAM4(4,1)=b0*4; yyirfAM4(:,1)=DD4+ZZ4*ssirfAM4(:,1);
for t=2:xx+1 %ZLB period
    ssirfAM4(:,t)=CC4(:,2)+TT4(:,:,2)*ssirfAM4(:,t-1);
    yyirfAM4(:,t)=DD4+ZZ4*ssirfAM4(:,t);
end
for t=xx+2:nirf %Coordination
    ssirfAM4(:,t)=CC4(:,1)+TT4(:,:,1)*ssirfAM4(:,t-1);
    yyirfAM4(:,t)=DD4+ZZ4*ssirfAM4(:,t);
end
%% Simulations Conditional On Regime Paths

ssirfshsp=zeros(size(TTshsp,1),nirf); yyirfshsp=zeros(size(ZZshsp,1),nirf);
ssirfshsp(4,1)=b0*4; ssirfshsp(12,1)=b0*4; yyirfshsp(:,1)=DDshsp+ZZshsp*ssirfshsp(:,1);
ssirf_highdebt=zeros(size(TTshsp,1),nirf); yyirf_highdebt=zeros(size(ZZshsp,1),nirf);
ssirf_highdebt(4,1)=(b0+.30)*4; ssirf_highdebt(12,1)=(b0+.30)*4;
yyirf_highdebt(:,1)=DDshsp+ZZshsp*ssirf_highdebt(:,1);

for t=2:xx+1 %negative demand shock period
    ssirfshsp(:,t)=CCshsp(:,2)+TTshsp(:,:)*ssirfshsp(:,t-1);
    yyirfshsp(:,t)=DDshsp+ZZshsp*ssirfshsp(:,t);
    ssirf_highdebt(:,t)=CCshsp(:,2)+TTshsp(:,:)*ssirf_highdebt(:,t-1);
    yyirf_highdebt(:,t)=DDshsp+ZZshsp*ssirf_highdebt(:,t);
end

for t=xx+2:nirf %end of recession
    ssirfshsp(:,t)=CCshsp(:,1)+TTshsp*ssirfshsp(:,t-1);
    yyirfshsp(:,t)=DDshsp+ZZshsp*ssirfshsp(:,t);
    ssirf_highdebt(:,t)=CCshsp(:,1)+TTshsp*ssirf_highdebt(:,t-1);
    yyirf_highdebt(:,t)=DDshsp+ZZshsp*ssirf_highdebt(:,t);
end
yyirfshsp(4,:)=yyirfshsp(4,:)/4;
yyirfshsp(12,:)=yyirfshsp(12,:)/4;
yyirfAM4(4,:)=yyirfAM4(4,:)/4;
yyirf_highdebt(4,:)=yyirf_highdebt(4,:)/4;

%% Computing Conditional Moments

nsb_shsp=length(HBss);
for ib=1:nsb_shsp
    TTshspH(:,:,ib)=[TTshsp(:,:),CCshsp(:,ib)];
    RQshspH(:,:,ib)=RRshsp(:,:)*QQshsp.^0.5;
end

[~,OMEGAHshsp,CSIHshsp] = mss_lom_nonhomog6_note_plus_dsge_RQ6_rev(TTshspH,RQshspH,HBss,0);

smbc_shsp=zeros(size(ssirfshsp,2)-1,nsb_shsp);
smbc_shsp(1:10,2)=ones(10,1); %LD (FA is expecting no conflict)
smbc_shsp(11:30,1)=ones(20,1); %Out of the LD (high demand)

if cond_moments==1
    [ ~, SD_shsp, ~, ~, ~, ~ ] = cond_moments93_co_pr_rev( OMEGAHshsp, CSIHshsp, nsb_shsp*1, Thor, 1, ssirfshsp(:,2:end)', smbc_shsp);
    save CondMoms\SDshsp SD_shsp
else
    load SDshsp SD_shsp
end
%%
posvar=[2;26]; % set the position of the variables in the quadratic function of welfare

[ welfareloss_shsp ] = welfare_msdsge_rev( OMEGAHshsp, CSIHshsp, nsb_shsp*1, 1, p_st, p_ss,oparams,posvar, ssirfshsp(:,2:end)', smbc_shsp );

% Welfare assuming lower beta 
p_sslb=p_ss; p_sslb(14)=oparams(6);
[ welfareloss_shsp_beta ] = welfare_msdsge_rev( OMEGAHshsp, CSIHshsp, nsb_shsp*1, 1, p_st, p_sslb,oparams,posvar, ssirfshsp(:,2:end)', smbc_shsp );


%% Producing Figures 6, 7, and 8

%% Figure 7
figure('name','Rankings of policies')
subplot(1,2,1)
plot(1:nirf-1,-welfarelossAF'*100,'--k')
hold on
plot(1:nirf-1,-welfarelossLCF'*100,'k')
hold on
plot(1:nirf-1,-welfarelossLCM'*100,'-ok','MarkerSize',3)
hold on
plot(1:nirf-1,-welfareloss_shsp'*100,'-pk')
hold on

if bw==1
    patch([1 xx xx 1],[min(min([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])])) ...
        min(min([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])]))...
        max(max([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])]))...
        max(max([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])])) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[min(min(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100]))) ...
        min(min(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])))...
        max(max(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])))...
        max(max(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])))],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
else
    patch([1 xx xx 1],[min(min([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])])) ...
        min(min([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])]))...
        max(max([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])]))...
        max(max([-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])])) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[min(min(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100]))) ...
        min(min(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])))...
        max(max(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])))...
        max(max(-([welfarelossAF'*100,welfarelossLCF'*100,welfarelossLCM'*100,welfareloss_shsp'*100])))],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
end
hold off
title('Baseline Beta')
axis tight

subplot(1,2,2)
plot(1:nirf-1,-welfarelossAFbeta'*100,'--k')
hold on
plot(1:nirf-1,-welfarelossLCFbeta'*100,'k')
hold on
plot(1:nirf-1,-welfarelossLCMbeta'*100,'o-k','MarkerSize',3)
hold on
plot(1:nirf-1,-welfareloss_shsp_beta'*100,'-pk')
hold on
if bw==1
    patch([1 xx xx 1],[min(min([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])])) ...
        min(min([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])]))...
        max(max([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])]))...
        max(max([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])])) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[min(min(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100]))) ...
        min(min(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])))...
        max(max(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])))...
        max(max(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])))],'k','EdgeColor','k','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
else
    patch([1 xx xx 1],[min(min([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])])) ...
        min(min([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])]))...
        max(max([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])]))...
        max(max([-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])])) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
    hold on
    patch([xx xx+ft xx+ft xx],[min(min(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100]))) ...
        min(min(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])))...
        max(max(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])))...
        max(max(-([welfarelossAFbeta'*100,welfarelossLCFbeta'*100,welfarelossLCMbeta'*100,welfareloss_shsp_beta'*100])))],'r','EdgeColor','r','FaceAlpha' ,.25,'EdgeAlpha' ,.25)
end

hold off
title('Low Beta')
h1=legend('LD->FL','LD->Conflict(FL)->FL','LD->Conflict(ML)->ML','Emergency Budget','Location','NorthWest');
axis tight
if dufig==1
    savefigure_pdf('Figures\Figure7') ;
end

%% Figure 8
hori=[4,20,40];
figure ('name','Evolution of Uncertainty (1y,5y,10y)')
for i=1:3
    
    subplot(3,2,1+2*(i-1))
    plot(1:nirf-1,100*squeeze(SD_AF(1,:,1+hori(i))),'--k')
    hold on
    plot(1:nirf-1,100*squeeze(SD_LCF(1,:,1+hori(i))),'k')
    hold on
    plot(1:nirf-1,100*squeeze(SD_LCM(1,:,1+hori(i))),'o-k','MarkerSize',3)
    hold on
    plot(1:nirf-1,100*squeeze(SD_shsp(1,:,1+hori(i))),'-pk')
    hold off
    axis tight
    ylim([5.5,8])
    yticks([5.5 6 6.5 7 7.5 8.0])
    if i==1
        h2=ylabel('One-Year Horizon');
        set(h2,'FontSize',12);
    elseif i==2
        h3=ylabel('Five-Year Horizon');
        set(h3,'FontSize',12);
    elseif i==3
        h4=ylabel('Ten-Year Horizon');
        set(h4,'FontSize',12);
    end
    ytickformat('%.2f')
    if i==1
        title('Output Gap Uncertainty')
    end
    if i==2
        h1=legend('LD->FL','LD->Conflict(FL)->FL','LD->Conflict(ML)->ML','Emergency Budget');
        set(h1,'FontSize',7)
    end
    
    subplot(3,2,2+2*(i-1))
    plot(1:nirf-1,100*squeeze(SD_AF(2,:,1+hori(i))),'--k')
    hold on
    plot(1:nirf-1,100*squeeze(SD_LCF(2,:,1+hori(i))),'k')
    hold on
    plot(1:nirf-1,100*squeeze(SD_LCM(2,:,1+hori(i))),'o-k','MarkerSize',3)
    hold on
    plot(1:nirf-1,100*squeeze(SD_shsp(2,:,1+hori(i))),'-pk')
    hold off
    axis tight
    ylim([0.35,0.6])
    yticks([.35 .4 .45 .5 .55 .6])
    if i==1
        h2=ylabel('One-Year Horizon');
        set(h2,'FontSize',12);
    elseif i==2
        h3=ylabel('Five-Year Horizon');
        set(h3,'FontSize',12);
    elseif i==3
        h4=ylabel('Ten-Year Horizon');
        set(h4,'FontSize',12);
    end
    ytickformat('%.2f')
    if i==1
        title('Inflation Uncertainty')
    end
end
if dufig==1
    savefigure_pdf('Figures\Figure8');
end

%% Figure 6
figure('name','Shock-Specific Rule Comparison with always AM (Benchmark) #2- Observables','numbertitle','off')

subplot(2,3,1)
if bw==1
    plot(0:nirf-1,ssirfshsp(1,:)*100,'k');
else
    plot(0:nirf-1,ssirfshsp(1,:)*100,'b');
end
hold on
patch([1 xx xx 1],[min([(ssirfshsp(1,:)*100),ssirfAM4(1,:)*100,0]) ...
    min([(ssirfshsp(1,:)*100),ssirfAM4(1,:)*100,0])...
    max([(ssirfshsp(1,:)*100),ssirfAM4(1,:)*100,0])...
    max([(ssirfshsp(1,:)*100),ssirfAM4(1,:)*100,0]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,ssirfAF(1,:)*0,'-.r','LineWidth',1);
hold off
axis tight
title('Output Gap')

subplot(2,3,2)
if bw==1
    plot(0:nirf-1,yyirfshsp(2,:)*400,'k');
else
    plot(0:nirf-1,yyirfshsp(2,:)*400,'b');
end
hold on
patch([1 xx xx 1],[min([yyirfshsp(2,:)*400,yyirfAM4(2,:)*40,1]) ...
    min([yyirfshsp(2,:)*400,yyirfAM4(2,:)*40,1])...
    max([yyirfshsp(2,:)*400,yyirfAM4(2,:)*40,1])...
    max([yyirfshsp(2,:)*400,yyirfAM4(2,:)*40,1]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,ones(nirf)*DD(2)*400,'-.r','LineWidth',1);
hold off
axis tight
title('Inflation')


subplot(2,3,3)
if bw==1
    plot(0:nirf-1,yyirfshsp(3,:)*400,'k');
else
    plot(0:nirf-1,yyirfshsp(3,:)*400,'b');
end
hold on
patch([1 xx xx 1],[min([yyirfshsp(3,:)*400,yyirfAM4(3,:)*400,0]) ...
    min([yyirfshsp(3,:)*400,yyirfAM4(3,:)*400,0])...
    max([yyirfshsp(3,:)*400,yyirfAM4(3,:)*400,0])...
    max([yyirfshsp(3,:)*400,yyirfAM4(3,:)*400,0]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
hold on
plot(0:nirf-1,ones(nirf)*DD(3)*400,'-.r','LineWidth',1);
hold off
axis tight
title('FFR')

subplot(2,3,4)
if bw==1
    plot(0:nirf-1,yyirfshsp(4,:)*100,'k');
    hold on
    plot(0:nirf-1,yyirfshsp(12,:)*100,'*--k','linewidth',1);
else
    plot(0:nirf-1,yyirfshsp(4,:)*100,'b');
    hold on
    plot(0:nirf-1,yyirfshsp(12,:)*100,'*--r','linewidth',1);
end
hold on

X=[0:nirf-1,fliplr(0:nirf-1)];                %#create continuous x value array for plotting
Y=[yyirfshsp(4,:)*100,fliplr(yyirfshsp(12,:)*100)];              %#create y values for out and then back
fill(X,Y,'b','EdgeColor','b','FaceAlpha' ,.25,'EdgeAlpha' ,.25);
hold on
patch([1 xx xx 1],[min([yyirfshsp(4,:)*100,yyirfshsp(12,:)*100,yyirfAM4(4,:)*100]) ...
    min([yyirfshsp(4,:)*100,yyirfshsp(12,:)*100,yyirfAM4(4,:)*100])...
    max([yyirfshsp(4,:)*100,yyirfshsp(12,:)*100,yyirfAM4(4,:)*100])...
    max([yyirfshsp(4,:)*100,yyirfshsp(12,:)*100,yyirfAM4(4,:)*100]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
hold off
axis tight
title('Debt-to-GDP Ratio')

h1=legend('Total Debt','Regular-Budget Debt','Emergency-Budget Debt');
set(h1,'Location','SouthWest','FontSize',10);

subplot(2,3,5)
if bw==1
    plot(0:nirf-1,yyirfshsp(9,:),'k');
else
    plot(0:nirf-1,yyirfshsp(9,:),'b');
end
hold on
patch([1 xx xx 1],[min([yyirfshsp(9,:),yyirfAM4(9,:)]) ...
    min([yyirfshsp(9,:),yyirfAM4(9,:)])...
    max([yyirfshsp(9,:),yyirfAM4(9,:)])...
    max([yyirfshsp(9,:),yyirfAM4(9,:)]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,ones(nirf)*DD(9),'-.r','LineWidth',1);
hold off
axis tight
title('Price of Bonds')


subplot(2,3,6)
if bw==1
    plot(0:nirf-1,yyirfshsp(8,:)*400,'k');
else
    plot(0:nirf-1,yyirfshsp(8,:)*400,'b');
end
hold on
patch([1 xx xx 1],[min([yyirfshsp(8,:)*400,yyirfAM4(8,:)*400,0]) ...
    min([yyirfshsp(8,:)*400,yyirfAM4(8,:)*400,0])...
    max([yyirfshsp(8,:)*400,yyirfAM4(8,:)*400,0])...
    max([yyirfshsp(8,:)*400,yyirfAM4(8,:)*400,0]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,ones(nirf)*DD(8)*400,'-.r','LineWidth',1);
hold off
axis tight
title('Real FFR')


if dufig==1
    savefigure_pdf('Figures\Figure6') ;
end


%% Bigger Initial Shock

p_st_scspLS=p_st_scsp;
p_st_scspLS(35,:)=p_st(35,[1,3]);
p_st_scspLS(35,2)=-.3;
b2=.07; 

HBbb=zeros(2,2); HBbb(1,1)=oparams(1); HBbb(2,1)=1-oparams(1); HBbb(1,2)=1-oparams(2); HBbb(2,2)=oparams(2);
[GAM0scspLS,GAM1scspLS,PSIscspLS,PPIscspLS,snames] = gammas_2771_shock_specific_Rules_4(p_st_scspLS,p_ss,msel,HBbb);
[TTshspLS,~,RRshspLS,fmat,fwt,ywt,gev,euLS,loose]=gensys(GAM0scspLS,GAM1scspLS,zeros(size(GAM0scspLS,1),1),PSIscspLS,PPIscspLS);

if sum(euLS)~=2
    error('indeterminacy when solving the model with shock-specific policy rules - Large Shock case')
end
%Remove the dummy variables
CCshspLS(:,1)=RRshspLS(:,10); CCshspLS(:,2)=RRshspLS(:,11);
%Resizing
CCshspLS=CCshspLS([1:23,26:size(CCshspLS,1)],:);
TTshspLS=TTshspLS([1:23,26:size(TTshspLS,1)],[1:23,26:size(TTshspLS,1)]);
RRshspLS=RRshspLS(:,1:9);
ZZshspLS=zeros(size(TTshspLS)); ZZshspLS(1,1:size(ZZ,2))=ZZ(1,:); ZZshspLS(2:16,2:16)=eye(15);
ZZshspLS(8,:)=ZZshspLS(8,:)*0; ZZshspLS(8,3)=1; ZZshspLS(8,27)=1;  %the real rate in the actual economy
ZZshspLS(9,:)=ZZshspLS(9,:)*0; ZZshspLS(9,7)=1;   %the price of long-term bonds

DDshspLS=[DD;zeros(size(TTshspLS,1)-length(DD),1)]; DDshspLS(10:13,1)=DD(1:4,1);DDshspLS(12)=DDshspLS(13);

HBlc5=HBlc4;
p_st_lc5=p_st_lc4;

%Bigger initial shock
p_st_lc5(35,2)=p_st_scspLS(35,end);
p_st_lc5(45,2)=(1-oparams(2))*p_st_lc5(35,1)+oparams(2)*p_st_lc5(35,2);
p_st_lc5(45,1)=oparams(1)*p_st_lc5(35,1)+(1-oparams(1))*p_st_lc5(35,2);


[TT5,QQ5,RR5,HH5,DD5,ZZ5,VV5,RC5,CC5,MATmo5]=...
    fb_sysmat3_bl8(p_st_lc5,p_ss,p_er,o_er,HBlc5,msel,MATlc);


if sum(RC5)~=1
    error('indeterminacy when solving the model with shock-specific policy rules - Active exit and big demand shock')
end


%impulse response
ssirfshspLS=zeros(size(TTshspLS,1),nirf); yyirfshspLS=zeros(size(ZZshspLS,1),nirf);
ssirfshspLS2=zeros(size(TTshspLS,1),nirf); yyirfshspLS2=zeros(size(ZZshspLS,1),nirf);
ssirf5=zeros(size(TT5,1),nirf); yyirf5=zeros(size(ZZ5,1),nirf);
% ssirfshsp(:,1)=CCshsp(:,1)+TTshsp(:,4)*b0*4+TTshsp(:,12)*b0*4;
ssirfshspLS(4,1)=b2*4; ssirfshspLS(12,1)=(b2)*4; yyirfshspLS(:,1)=DDshspLS+ZZshspLS*ssirfshspLS(:,1);
ssirfshspLS2(4,1)=b2*4; ssirfshspLS2(12,1)=(b2)*4; yyirfshspLS2(:,1)=DDshspLS+ZZshspLS*ssirfshspLS2(:,1);
ssirf5(4,1)=b2*4; yyirf5(:,1)=DD5+ZZ5*ssirf5(:,1);

for t=2:xx+1 %negative demand shock period
    ssirfshspLS(:,t)=CCshspLS(:,2)+TTshspLS(:,:)*ssirfshspLS(:,t-1);
    yyirfshspLS(:,t)=DDshspLS+ZZshspLS*ssirfshspLS(:,t);
    ssirfshspLS2(:,t)=CCshspLS(:,2)+TTshspLS(:,:)*ssirfshspLS2(:,t-1);
    yyirfshspLS2(:,t)=DDshspLS+ZZshspLS*ssirfshspLS2(:,t);
    ssirf5(:,t)=CC5(:,2)+TT5(:,:,2)*ssirf5(:,t-1);
    yyirf5(:,t)=DD5+ZZ5*ssirf5(:,t);
end

for t=xx+2:nirf %end of recession
    ssirfshspLS(:,t)=CCshspLS(:,1)+TTshspLS*ssirfshspLS(:,t-1);
    yyirfshspLS(:,t)=DDshspLS+ZZshspLS*ssirfshspLS(:,t);
    ssirfshspLS2(:,t)=CCshspLS(:,1)+TTshspLS*ssirfshspLS2(:,t-1);
    yyirfshspLS2(:,t)=DDshspLS+ZZshspLS*ssirfshspLS2(:,t);
    ssirf5(:,t)=CC5(:,1)+TT5(:,:,1)*ssirf5(:,t-1);
    yyirf5(:,t)=DD5+ZZ5*ssirf5(:,t);
end

%rescaling the debt-to-GDP ratio
yyirfshspLS(4,:)=yyirfshspLS(4,:)/4;
yyirfshspLS2(4,:)=yyirfshspLS2(4,:)/4;
yyirf5(4,:)=yyirf5(4,:)/4;
%% Figure 9
figure('name','Shock-Specific Rule: Large Demand Shock - Observables','numbertitle','off')

subplot(2,3,1)
if bw==1
    plot(0:nirf-1,ssirfshspLS(1,:)*100,'k');
else
    plot(0:nirf-1,ssirfshspLS(1,:)*100,'b');
end
hold on
plot(0:nirf-1,ssirf5(1,:)*100,'--k');
hold on
patch([1 xx xx 1],[min([(ssirfshspLS(1,:)*100),ssirf5(1,:)*100,0]) ...
    min([(ssirfshspLS(1,:)*100),ssirf5(1,:)*100,0])...
    max([ssirfshspLS(1,:)*100,ssirf5(1,:)*100,0])...
    max([ssirfshspLS(1,:)*100,ssirf5(1,:)*100,0]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,yyirfshspLS*0,'-.r','LineWidth',1);
hold off
axis tight
title('Output Gap')


subplot(2,3,2)
if bw==1
    plot(0:nirf-1,yyirfshspLS(2,:)*400,'k');
else
    plot(0:nirf-1,yyirfshspLS(2,:)*400,'b');
end
hold on
plot(0:nirf-1,yyirf5(2,:)*400,'--k');
hold on
patch([1 xx xx 1],[min([(yyirfshspLS(2,:)*400),yyirf5(2,:)*400,DDshspLS(2)*400,1]) ...
    min([(yyirfshspLS(2,:)*400),yyirf5(2,:)*400,DDshspLS(2)*400,1])...
    max([yyirfshspLS(2,:)*400,yyirf5(2,:)*400,DDshspLS(2)*400,1])...
    max([yyirfshspLS(2,:)*400,yyirf5(2,:)*400,DDshspLS(2)*400,1]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
hold on
plot(0:nirf-1,DDshspLS(2)*400*ones(nirf,1),'-.r','LineWidth',1);
hold off
axis tight
title('Inflation')


subplot(2,3,3)
if bw==1
    plot(0:nirf-1,yyirfshspLS(3,:)*400,'k');
else
    plot(0:nirf-1,yyirfshspLS(3,:)*400,'b');
end
hold on
plot(0:nirf-1,yyirf5(3,:)*400,'--k');
hold on
patch([1 xx xx 1],[min([(yyirfshspLS(3,:)*400),yyirf5(3,:)*400,0]) ...
    min([(yyirfshspLS(3,:)*400),yyirf5(3,:)*400,0])...
    max([yyirfshspLS(3,:)*400,yyirf5(3,:)*400,0])...
    max([yyirfshspLS(3,:)*400,yyirf5(3,:)*400,0]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,yyirfshspLS(1,:)*0,'-.r','LineWidth',1);
hold off
axis tight
title('FFR')


subplot(2,3,4)
if bw==1
    plot(0:nirf-1,yyirfshspLS(4,:)*100,'k');
else
    plot(0:nirf-1,yyirfshspLS(4,:)*100,'b');
end
hold on
plot(0:nirf-1,yyirf5(4,:)*100,'--k');
hold on
patch([1 xx xx 1],[min([(yyirfshspLS(4,:)*100),(yyirf5(4,:)*100),DDshspLS(4)*25]) ...
    min([(yyirfshspLS(4,:)*100),(yyirf5(4,:)*100),DDshspLS(4)*25])...
    max([yyirfshspLS(4,:)*100,(yyirf5(4,:)*100),DDshspLS(4)*25])...
    max([yyirfshspLS(4,:)*100,(yyirf5(4,:)*100),DDshspLS(4)*25]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,DDshspLS(4)*25*ones(nirf,1),'-.r','LineWidth',1);
hold off
axis tight
title('Debt-to-GDP Ratio')
h1=legend('Emergency-Budget Rule','Always ML Rule');
set(h1,'Location','SouthWest','FontSize',10);

subplot(2,3,5)
if bw==1
    plot(0:nirf-1,yyirfshspLS(9,:),'k');
else
    plot(0:nirf-1,yyirfshspLS(9,:),'b');
end
hold on
plot(0:nirf-1,yyirf5(9,:),'--k');
hold on
patch([1 xx xx 1],[min([(yyirfshspLS(9,:)),yyirf5(9,:)]) ...
    min([(yyirfshspLS(9,:)),yyirf5(9,:)])...
    max([yyirfshspLS(9,:),yyirf5(9,:)])...
    max([yyirfshspLS(9,:),yyirf5(9,:)]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,DDshspLS(9)*ones(nirf,1),'-.r','LineWidth',1);
hold off
axis tight
title('Price of Bonds')


subplot(2,3,6)
if bw==1
    plot(0:nirf-1,yyirfshspLS(8,:)*400,'k');
else
    plot(0:nirf-1,yyirfshspLS(8,:)*400,'b');
end
hold on
plot(0:nirf-1,yyirf5(8,:)*400,'--k');
hold on
patch([1 xx xx 1],[min([(yyirfshspLS(8,:)*400),yyirf5(8,:)*400,0]) ...
    min([(yyirfshspLS(8,:)*400),yyirf5(8,:)*400,0])...
    max([yyirfshspLS(8,:)*400,yyirf5(8,:)*400,0])...
    max([yyirfshspLS(8,:)*400,yyirf5(8,:)*400,0]) ],'k','EdgeColor','k','FaceAlpha' ,.5,'EdgeAlpha' ,.5)
hold on
plot(0:nirf-1,yyirfshspLS(1,:)*0,'-.r','LineWidth',1);
hold off
axis tight
title('Real FFR')

if dufig==1
    savefigure_pdf('Figures\Figure9') ;
end
