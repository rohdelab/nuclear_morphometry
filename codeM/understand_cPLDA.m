clc
clear all
close all
CL={[0.5 0.5 0.5],[0.5 0 0.5],[0.75 0.75 0],'m',[0.4 0.4 0.4],[0.93 0.45 0.10],[0.27 0.62 0.25],[0.79 0.15 0.18],[0.09 0.44 0.65]};

%% PARAMETERS
alpha=1.1618; % 1.1618; 0.00

%% INPUT DATASETS
X1 = [0.5893    0.6515    0.7851    0.9165    0.8635    0.6262 0.3934    0.3427    0.4510    0.5962    0.5870    0.7921 0.9487    0.6561    0.4165    0.0939    0.0616    0.1884 0.2713    0.0432    0.1999    0.4441    0.4879    0.2990 0.2252    0.1469    0.0962    0.4234    0.4303    0.2368
    0.8153    0.8883    0.8504    0.7161    0.5847    0.5672 0.6431    0.8416    0.9234    0.7803    0.9788    0.6839 0.4066    0.5321    0.5993    0.1000    0.2985    0.2898 0.4504    0.5117    0.2693    0.0825    0.3394    0.4328 0.2139    0.4328    0.6752    0.3073    0.4737    0.5088 ];
y1=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

X2=X1; X2(1,:)=-X2(1,:);
y2=y1;

% % % uncomment originally
% anggg=20;
% X1=[cosd(anggg) -sind(anggg);sind(anggg) cosd(anggg)]*[0.75 0;0 0.4]*X1;
% X2=[cosd(anggg) -sind(anggg);sind(anggg) cosd(anggg)]*[0.25 0;0 0.5]*X2;

X1_m0=X1-repmat(mean(X1')',1,size(X1,2));
X2_m0=X2-repmat(mean(X2')',1,size(X2,2));

Xall=[X1 X2];
Alllabel=[y1 y2];

Xall_m0=Xall-repmat(mean(Xall')',1,size(Xall,2));

xcvxcv
%%
% Feature scaling
clear ss1;
for aa=1:size(X1_m0,1)
    ss1(aa,aa)=1/norm(X1_m0(aa,:));
end
% X1_m0=ss1*X1_m0;

clear ss2;
for aa=1:size(X2_m0,1)
    ss2(aa,aa)=1/norm(X2_m0(aa,:));
end
% X2_m0=ss2*X2_m0;
%

% Feature scaling
clear ssall;
for aa=1:size(Xall_m0,1)
    ssall(aa,aa)=1/norm(Xall_m0(aa,:));
end
% Xall_m0=ssall*Xall_m0;
%


%%
[bpca_org_1,~,~]=svd(X1_m0,'econ');
[bplda_pca_1,~] = PLDA(bpca_org_1'*X1_m0,y1,alpha); bplda_pca_1=real(bplda_pca_1);
bplda_org_1=bpca_org_1*bplda_pca_1;

bplda_org_1(:,1)'*bplda_org_1(:,2)

[bpca_org_2,~,~]=svd(X2_m0,'econ');
[bplda_pca_2,~] = PLDA(bpca_org_2'*X2_m0,y2,alpha); bplda_pca_2=real(bplda_pca_2);
bplda_org_2=bpca_org_2*bplda_pca_2; bplda_org_2(:,1)=-bplda_org_2(:,1);

bplda_org_2(:,1)'*bplda_org_2(:,2)
%

[bpca_org_all,~,~]=svd(Xall_m0,'econ');
[ST1,SWNew1,ss1] = cPLDA_p1(bpca_org_all'*X1_m0,y1,alpha);
[ST2,SWNew2,ss2] = cPLDA_p1(bpca_org_all'*X2_m0,y2,alpha);


ST=ST1+ST2; SWNew=SWNew1+SWNew2;

Ndir=2;
[bplda_pca_all,~] = cPLDA_p2(ST,SWNew,Ndir,ss1); % or ss2 or a garbage
bplda_org_all=bpca_org_all*bplda_pca_all;

bplda_org_all(:,1)'*bplda_org_all(:,2)

%% A
X_org_1=X1_m0+ss1*repmat(mean(X1')',1,size(X1,2));
X_org_2=X2_m0+ss2*repmat(mean(X2')',1,size(X2,2));

dir1plda_org_1=bplda_org_1(:,1); dir1plda_org_1=(1)*dir1plda_org_1;
dir2plda_org_1=bplda_org_1(:,2); dir2plda_org_1=(0.5)*dir2plda_org_1;

dir1plda_org_2=bplda_org_2(:,1); dir1plda_org_2=(1)*dir1plda_org_2;
dir2plda_org_2=bplda_org_2(:,2); dir2plda_org_2=(0.5)*dir2plda_org_2;

dir1plda_org_all=bplda_org_all(:,1); dir1plda_org_all=(1)*dir1plda_org_all;
dir2plda_org_all=bplda_org_all(:,2); dir2plda_org_all=(0.5)*dir2plda_org_all;

%% B
X_org_1_b=inv(ss1)*X1_m0+(1)*repmat(mean(X1')',1,size(X1,2));
X_org_2_b=inv(ss2)*X2_m0+(1)*repmat(mean(X2')',1,size(X2,2));

bplda_org_1_b=inv(ss1)*bplda_org_1;
bplda_org_2_b=inv(ss1)*bplda_org_2;
bplda_org_all_b=inv(ss1)*bplda_org_all;
for i=1:size(bplda_org_1_b,2)
    bplda_org_1_b(:,i)=bplda_org_1_b(:,i)/sqrt(bplda_org_1_b(:,i)'*bplda_org_1_b(:,i));
end
for i=1:size(bplda_org_2_b,2)
    bplda_org_2_b(:,i)=bplda_org_2_b(:,i)/sqrt(bplda_org_2_b(:,i)'*bplda_org_2_b(:,i));
end
for i=1:size(bplda_org_all_b,2)
    bplda_org_all_b(:,i)=bplda_org_all_b(:,i)/sqrt(bplda_org_all_b(:,i)'*bplda_org_all_b(:,i));
end

dir1plda_org_1_b=bplda_org_1_b(:,1); dir1plda_org_1_b=(1)*dir1plda_org_1_b;
dir2plda_org_1_b=bplda_org_1_b(:,2); dir2plda_org_1_b=(0.5)*dir2plda_org_1_b;

dir1plda_org_2_b=bplda_org_2_b(:,1); dir1plda_org_2_b=(1)*dir1plda_org_2_b;
dir2plda_org_2_b=bplda_org_2_b(:,2); dir2plda_org_2_b=(0.5)*dir2plda_org_2_b;

dir1plda_org_all_b=bplda_org_all_b(:,1); dir1plda_org_all_b=(1)*dir1plda_org_all_b;
dir2plda_org_all_b=bplda_org_all_b(:,2); dir2plda_org_all_b=(0.5)*dir2plda_org_all_b;

%% C
X_org_1_c=inv(ss1)*X1_m0+(1)*repmat(mean(X1')',1,size(X1,2));
X_org_2_c=inv(ss2)*X2_m0+(1)*repmat(mean(X2')',1,size(X2,2));

bplda_org_1_c=inv(ss2)*bplda_org_1;
bplda_org_2_c=inv(ss2)*bplda_org_2;
bplda_org_all_c=inv(ss2)*bplda_org_all;
for i=1:size(bplda_org_1_c,2)
    bplda_org_1_c(:,i)=bplda_org_1_c(:,i)/sqrt(bplda_org_1_c(:,i)'*bplda_org_1_c(:,i));
end
for i=1:size(bplda_org_2_c,2)
    bplda_org_2_c(:,i)=bplda_org_2_c(:,i)/sqrt(bplda_org_2_c(:,i)'*bplda_org_2_c(:,i));
end
for i=1:size(bplda_org_all_c,2)
    bplda_org_all_c(:,i)=bplda_org_all_c(:,i)/sqrt(bplda_org_all_c(:,i)'*bplda_org_all_c(:,i));
end

dir1plda_org_1_c=bplda_org_1_c(:,1); dir1plda_org_1_c=(1)*dir1plda_org_1_c;
dir2plda_org_1_c=bplda_org_1_c(:,2); dir2plda_org_1_c=(0.5)*dir2plda_org_1_c;

dir1plda_org_2_c=bplda_org_2_c(:,1); dir1plda_org_2_c=(1)*dir1plda_org_2_c;
dir2plda_org_2_c=bplda_org_2_c(:,2); dir2plda_org_2_c=(0.5)*dir2plda_org_2_c;

dir1plda_org_all_c=bplda_org_all_c(:,1); dir1plda_org_all_c=(1)*dir1plda_org_all_c;
dir2plda_org_all_c=bplda_org_all_c(:,2); dir2plda_org_all_c=(0.5)*dir2plda_org_all_c;

%% D
X_org_1_d=inv(ss1)*X1_m0+(1)*repmat(mean(X1')',1,size(X1,2));
X_org_2_d=inv(ss2)*X2_m0+(1)*repmat(mean(X2')',1,size(X2,2));

bplda_org_1_d=inv(ssall)*bplda_org_1;
bplda_org_2_d=inv(ssall)*bplda_org_2;
bplda_org_all_d=inv(ssall)*bplda_org_all;
for i=1:size(bplda_org_1_d,2)
    bplda_org_1_d(:,i)=bplda_org_1_d(:,i)/sqrt(bplda_org_1_d(:,i)'*bplda_org_1_d(:,i));
end
for i=1:size(bplda_org_2_d,2)
    bplda_org_2_d(:,i)=bplda_org_2_d(:,i)/sqrt(bplda_org_2_d(:,i)'*bplda_org_2_d(:,i));
end
for i=1:size(bplda_org_all_d,2)
    bplda_org_all_d(:,i)=bplda_org_all_d(:,i)/sqrt(bplda_org_all_d(:,i)'*bplda_org_all_d(:,i));
end


dir1plda_org_1_d=bplda_org_1_d(:,1); dir1plda_org_1_d=(1)*dir1plda_org_1_d;
dir2plda_org_1_d=bplda_org_1_d(:,2); dir2plda_org_1_d=(0.5)*dir2plda_org_1_d;

dir1plda_org_2_d=bplda_org_2_d(:,1); dir1plda_org_2_d=(1)*dir1plda_org_2_d;
dir2plda_org_2_d=bplda_org_2_d(:,2); dir2plda_org_2_d=(0.5)*dir2plda_org_2_d;

dir1plda_org_all_d=bplda_org_all_d(:,1); dir1plda_org_all_d=(1)*dir1plda_org_all_d;
dir2plda_org_all_d=bplda_org_all_d(:,2); dir2plda_org_all_d=(0.5)*dir2plda_org_all_d;

%%
f1=figure('position',[0 0 1000 800],'name',['Linear analysis on data']); movegui(f1,'northeast'); 
subplot 221; hold on; grid on
plot([0 0],[-1.2 1.2],'k--'); plot([-1.2 1.2],[0 0],'k--');

plot(X_org_1(1,find(y1==1)),X_org_1(2,find(y1==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(X_org_1(1,find(y1==2)),X_org_1(2,find(y1==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dir1plda_org_1(1),dir1plda_org_1(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_1(1),dir2plda_org_1(2),'--','color',CL{3},'linewidth',2);

plot(X_org_2(1,find(y2==1)),X_org_2(2,find(y2==1)),'v','markerfacecolor',CL{8},'markeredgecolor',CL{8})
plot(X_org_2(1,find(y2==2)),X_org_2(2,find(y2==2)),'o','markerfacecolor',CL{9},'markeredgecolor',CL{9})
h(2)=quiver(0,0,dir1plda_org_2(1),dir1plda_org_2(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_2(1),dir2plda_org_2(2),'--','color',CL{4},'linewidth',2);

h(3)=quiver(0,0,dir1plda_org_all(1),dir1plda_org_all(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_all(1),dir2plda_org_all(2),'--','color',CL{2},'linewidth',2,'maxheadsize',.4);

legend(h(1:3),'plda - 1','plda - 2','c-plda - all (1 & 2)','location','southwest')
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('Original space - feature scaled','interpreter','latex')
set(gca,'fontsize',14); axis square

subplot 222; hold on; grid on
plot([0 0],[-1.2 1.2],'k--'); plot([-1.2 1.2],[0 0],'k--');

plot(X_org_1_b(1,find(y1==1)),X_org_1_b(2,find(y1==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(X_org_1_b(1,find(y1==2)),X_org_1_b(2,find(y1==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dir1plda_org_1_b(1),dir1plda_org_1_b(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_1_b(1),dir2plda_org_1_b(2),'--','color',CL{3},'linewidth',2);

plot(X_org_2_b(1,find(y2==1)),X_org_2_b(2,find(y2==1)),'v','markerfacecolor',CL{8},'markeredgecolor',CL{8})
plot(X_org_2_b(1,find(y2==2)),X_org_2_b(2,find(y2==2)),'o','markerfacecolor',CL{9},'markeredgecolor',CL{9})
h(2)=quiver(0,0,dir1plda_org_2_b(1),dir1plda_org_2_b(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_2_b(1),dir2plda_org_2_b(2),'--','color',CL{4},'linewidth',2);

h(3)=quiver(0,0,dir1plda_org_all_b(1),dir1plda_org_all_b(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_all_b(1),dir2plda_org_all_b(2),'--','color',CL{2},'linewidth',2,'maxheadsize',.4);

xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('To dataset - 1','interpreter','latex')
set(gca,'fontsize',14); axis square


subplot 223; hold on; grid on
plot([0 0],[-1.2 1.2],'k--'); plot([-1.2 1.2],[0 0],'k--');

plot(X_org_1_c(1,find(y1==1)),X_org_1_c(2,find(y1==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(X_org_1_c(1,find(y1==2)),X_org_1_c(2,find(y1==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dir1plda_org_1_c(1),dir1plda_org_1_c(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_1_c(1),dir2plda_org_1_c(2),'--','color',CL{3},'linewidth',2);

plot(X_org_2_c(1,find(y2==1)),X_org_2_c(2,find(y2==1)),'v','markerfacecolor',CL{8},'markeredgecolor',CL{8})
plot(X_org_2_c(1,find(y2==2)),X_org_2_c(2,find(y2==2)),'o','markerfacecolor',CL{9},'markeredgecolor',CL{9})
h(2)=quiver(0,0,dir1plda_org_2_c(1),dir1plda_org_2_c(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_2_c(1),dir2plda_org_2_c(2),'--','color',CL{4},'linewidth',2);

h(3)=quiver(0,0,dir1plda_org_all_c(1),dir1plda_org_all_c(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_all_c(1),dir2plda_org_all_c(2),'--','color',CL{2},'linewidth',2,'maxheadsize',.4);

xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('To dataset - 2','interpreter','latex')
set(gca,'fontsize',14); axis square


subplot 224; hold on; grid on
plot([0 0],[-1.2 1.2],'k--'); plot([-1.2 1.2],[0 0],'k--');

plot(X_org_1_d(1,find(y1==1)),X_org_1_d(2,find(y1==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(X_org_1_d(1,find(y1==2)),X_org_1_d(2,find(y1==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dir1plda_org_1_d(1),dir1plda_org_1_d(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_1_d(1),dir2plda_org_1_d(2),'--','color',CL{3},'linewidth',2);

plot(X_org_2_d(1,find(y2==1)),X_org_2_d(2,find(y2==1)),'v','markerfacecolor',CL{8},'markeredgecolor',CL{8})
plot(X_org_2_d(1,find(y2==2)),X_org_2_d(2,find(y2==2)),'o','markerfacecolor',CL{9},'markeredgecolor',CL{9})
h(2)=quiver(0,0,dir1plda_org_2_d(1),dir1plda_org_2_d(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_2_d(1),dir2plda_org_2_d(2),'--','color',CL{4},'linewidth',2);

h(3)=quiver(0,0,dir1plda_org_all_d(1),dir1plda_org_all_d(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dir2plda_org_all_d(1),dir2plda_org_all_d(2),'--','color',CL{2},'linewidth',2,'maxheadsize',.4);

xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('To dataset - all','interpreter','latex')
set(gca,'fontsize',14); axis square




