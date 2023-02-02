clc
clear all
close all
CL={[0.5 0.5 0.5],[0.5 0 0.5],[0.75 0.75 0],'m',[0.4 0.4 0.4],[0.93 0.45 0.10],[0.27 0.62 0.25],[0.79 0.15 0.18],[0.09 0.44 0.65]};

%% PARAMETERS
alpha=1.618; %1.618; % 1.618; 0.00

%% INPUT DATASETS
X1 = [0.5893    0.6515    0.7851    0.9165    0.8635    0.6262 0.3934    0.3427    0.4510    0.5962    0.5870    0.7921 0.9487    0.6561    0.4165    0.0939    0.0616    0.1884 0.2713    0.0432    0.1999    0.4441    0.4879    0.2990 0.2252    0.1469    0.0962    0.4234    0.4303    0.2368
    0.8153    0.8883    0.8504    0.7161    0.5847    0.5672 0.6431    0.8416    0.9234    0.7803    0.9788    0.6839 0.4066    0.5321    0.5993    0.1000    0.2985    0.2898 0.4504    0.5117    0.2693    0.0825    0.3394    0.4328 0.2139    0.4328    0.6752    0.3073    0.4737    0.5088 ];
y1=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

X2=X1; X2(1,:)=-X2(1,:);
y2=y1; % fliplr(y1);


% % % uncomment originally
anggg=20;
X1=[cosd(anggg) -sind(anggg);sind(anggg) cosd(anggg)]*[0.75 0;0 0.4]*X1;
X2=[cosd(anggg) -sind(anggg);sind(anggg) cosd(anggg)]*[0.25 0;0 0.5]*X2;

xxT=[X1 X2];
lT=[y1 y2];

%%
uB{1}=X1; uB{2}=X2;
labelB{1}=y1; labelB{2}=y2;
%%
for a=1:length(uB)
    u=uB{a};
    [M2,K]=size(u);
    temp=u;
    PsiB{a}=temp;
    Psimean=mean(temp')';
    temp=temp-repmat(Psimean,1,K);
    PsiB_m0{a}=temp;
end


%%
AllData=xxT; lN=lT;
AllData_m0=AllData-repmat(mean(AllData')',1,size(AllData,2));
[Uall,~,~]=svd(AllData_m0,'econ');
clear ssall; temp=Uall'*AllData_m0;
for aa=1:size(temp,1)
    ssall(aa,aa)=1/norm(temp(aa,:));
end
for a=1:length(PsiB_m0)
    Psi=PsiB_m0{a};
    PsiPCA=Uall'*Psi;
    PsiPCAB{a}=PsiPCA;
end
[bplda_pca_joint,~] = cPLDA(PsiPCAB,labelB,alpha,ssall); bplda_pca_joint=real(bplda_pca_joint);
bplda_org_joint=Uall*bplda_pca_joint;
bplda_org_joint(:,1)'*bplda_org_joint(:,2)

[bplda_pca_joint2,~] = cPLDA_v2(PsiPCAB,labelB,alpha); bplda_pca_joint2=real(bplda_pca_joint2);
bplda_org_joint2=Uall*bplda_pca_joint2;
bplda_org_joint2(:,1)'*bplda_org_joint2(:,2)

[bplda_pca_1,~] = PLDA(PsiPCAB{1},labelB{1},alpha); bplda_pca_1=real(bplda_pca_1);
bplda_org_1=Uall*bplda_pca_1;
bplda_org_1(:,1)'*bplda_org_1(:,2)

[bplda_pca_2,~] = PLDA(PsiPCAB{2},labelB{2},alpha); bplda_pca_2=real(bplda_pca_2);
bplda_org_2=Uall*bplda_pca_2;
bplda_org_2(:,1)'*bplda_org_2(:,2)


X_org_1=X1;
X_org_2=X2;

%%
dir1plda_org_1=bplda_org_1(:,1); dir1plda_org_1=(1)*dir1plda_org_1;
dir2plda_org_1=bplda_org_1(:,2); dir2plda_org_1=(0.5)*dir2plda_org_1;

dir1plda_org_2=bplda_org_2(:,1); dir1plda_org_2=(1)*dir1plda_org_2;
dir2plda_org_2=bplda_org_2(:,2); dir2plda_org_2=(0.5)*dir2plda_org_2;

dir1plda_org_joint=bplda_org_joint(:,1); dir1plda_org_joint=(1)*dir1plda_org_joint;
dir2plda_org_joint=bplda_org_joint(:,2); dir2plda_org_joint=(0.5)*dir2plda_org_joint;

dir1plda_org_joint2=bplda_org_joint2(:,1); dir1plda_org_joint2=(1)*dir1plda_org_joint2;
dir2plda_org_joint2=bplda_org_joint2(:,2); dir2plda_org_joint2=(0.5)*dir2plda_org_joint2;

%%
f1=figure('position',[0 0 1000 800],'name',['Linear analysis on data']); movegui(f1,'northeast');
hold on; grid on
plot([0 0],[-1.2 1.2],'k--'); plot([-1.2 1.2],[0 0],'k--');

plot(X_org_1(1,find(y1==1)),X_org_1(2,find(y1==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(X_org_1(1,find(y1==2)),X_org_1(2,find(y1==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dir1plda_org_1(1),dir1plda_org_1(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_1(1),dir2plda_org_1(2),'-','color',CL{3},'linewidth',2);

plot(X_org_2(1,find(y2==1)),X_org_2(2,find(y2==1)),'o','markerfacecolor',CL{8},'markeredgecolor',CL{8})
plot(X_org_2(1,find(y2==2)),X_org_2(2,find(y2==2)),'v','markerfacecolor',CL{9},'markeredgecolor',CL{9})
h(2)=quiver(0,0,dir1plda_org_2(1),dir1plda_org_2(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_2(1),dir2plda_org_2(2),'-','color',CL{4},'linewidth',2);

h(3)=quiver(0,0,dir1plda_org_joint(1),dir1plda_org_joint(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_joint(1),dir2plda_org_joint(2),'-','color',CL{2},'linewidth',2,'maxheadsize',.4);

h(4)=quiver(0,0,dir1plda_org_joint2(1),dir1plda_org_joint2(2),'color',CL{1},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_joint2(1),dir2plda_org_joint2(2),'-','color',CL{1},'linewidth',2,'maxheadsize',.4);


legend(h(1:4),'plda - 1','plda - 2','c-plda - v1 (1 & 2)','c-plda - v2 (1 & 2)','location','southwest')
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('Original space - feature scaled','interpreter','latex')
set(gca,'fontsize',14); %axis square

