clc
clear all
close all

%%
tag=1;

alpha=0; % 1.1618;
do_spca=1; % 1 yes do scpa, 0 no do generic pca in scpa style

%%
switch tag
    case 1
        Data = [ 0.5893    0.6515    0.7851    0.9165    0.8635    0.6262 0.3934    0.3427    0.4510    0.5962    0.5870    0.7921 0.9487    0.6561    0.4165    0.0939    0.0616    0.1884 0.2713    0.0432    0.1999    0.4441    0.4879    0.2990 0.2252    0.1469    0.0962    0.4234    0.4303    0.2368
            0.8153    0.8883    0.8504    0.7161    0.5847    0.5672 0.6431    0.8416    0.9234    0.7803    0.9788    0.6839 0.4066    0.5321    0.5993    0.1000    0.2985    0.2898 0.4504    0.5117    0.2693    0.0825    0.3394    0.4328 0.2139    0.4328    0.6752    0.3073    0.4737    0.5088 ];
    case 2
        Data = [ 0.3082    0.3404    0.4050    0.4856    0.5386    0.6146    0.6722    0.7298    0.7944    0.8888    0.9234    0.8335    0.6377    0.5593    0.4372    0.0985    0.1469    0.1999    0.2391    0.2851    0.3220    0.3543 0.3911    0.4234    0.4695    0.4971    0.5225    0.5570    0.5893    0.3934
            0.9321    0.8912    0.8270    0.7715    0.7044    0.6489    0.5847    0.5263    0.4504    0.3949    0.3248    0.4036    0.5409    0.6927    0.7540    0.8445    0.7628    0.7131    0.6547    0.5788    0.5263    0.4650 0.4328    0.3745    0.3336    0.2635    0.2139    0.1788    0.1146    0.3482 ];
    otherwise
        Data=rand(2,30);
end
label=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
%%
Data=Data*1;
Data0=mean(Data')';
Data=Data-repmat(Data0,1,size(Data,2));

%%
CL={[0.5 0.5 0.5],[0.5 0 0.5],[0.75 0.75 0],'m',[0.4 0.4 0.4],[0.93 0.45 0.10],[0.27 0.62 0.25],[0.79 0.15 0.18],[0.09 0.44 0.65]};
[basisPCA_org,~,~]=svd(Data,'econ'); basisPCA_org(:,1)=-basisPCA_org(:,1);

[VecPLDA_pca,~] = PLDA(basisPCA_org'*Data,label,alpha); VecPLDA_pca=real(VecPLDA_pca); VecPLDA_pca(:,1)=-VecPLDA_pca(:,1);
basisPLDA_org=basisPCA_org*VecPLDA_pca;

c1=1;c2=2;
Data10=Data(:,label==c1); Data20=Data(:,label==c2);
if do_spca==1
    Data10=Data10-repmat(mean(Data(:,label==c1)')',1,length(find(label==c1)));
    Data20=Data20-repmat(mean(Data(:,label==c2)')',1,length(find(label==c2)));
else
    Data10=Data10-repmat(mean(Data(:,label==c1|label==c2)')',1,length(find(label==c1)));
    Data20=Data20-repmat(mean(Data(:,label==c1|label==c2)')',1,length(find(label==c2)));
end

cm=Data10*Data10'+Data20*Data20';
[VecSPCA_org,~]=eig(cm); basisSPCA_org=VecSPCA_org(:,[2 1]);

%%
Data_org=Data;
dirORG1_org=[1;0]; dirORG1_org=(1)*dirORG1_org;
dirORG2_org=[0;1]; dirORG2_org=(0.5)*dirORG2_org;
dirPCA1_org=basisPCA_org(:,1); dirPCA1_org=(1)*dirPCA1_org;
dirPCA2_org=basisPCA_org(:,2); dirPCA2_org=(0.5)*dirPCA2_org;
dirPLDA1_org=basisPLDA_org(:,1); dirPLDA1_org=(1)*dirPLDA1_org;
dirPLDA2_org=basisPLDA_org(:,2); dirPLDA2_org=(0.5)*dirPLDA2_org;
dirSPCA1_org=basisSPCA_org(:,1); dirSPCA1_org=(1)*dirSPCA1_org;
dirSPCA2_org=basisSPCA_org(:,2); dirSPCA2_org=(0.5)*dirSPCA2_org;

Data_pca=basisPCA_org'*Data_org;
dirORG1_pca=basisPCA_org'*dirORG1_org;
dirORG2_pca=basisPCA_org'*dirORG2_org;
dirPCA1_pca=basisPCA_org'*dirPCA1_org;
dirPCA2_pca=basisPCA_org'*dirPCA2_org;
dirPLDA1_pca=basisPCA_org'*dirPLDA1_org;
dirPLDA2_pca=basisPCA_org'*dirPLDA2_org;
dirSPCA1_pca=basisPCA_org'*dirSPCA1_org;
dirSPCA2_pca=basisPCA_org'*dirSPCA2_org;

Data_plda=basisPLDA_org'*Data_org;
dirORG1_plda=basisPLDA_org'*dirORG1_org;
dirORG2_plda=basisPLDA_org'*dirORG2_org;
dirPCA1_plda=basisPLDA_org'*dirPCA1_org;
dirPCA2_plda=basisPLDA_org'*dirPCA2_org;
dirPLDA1_plda=basisPLDA_org'*dirPLDA1_org;
dirPLDA2_plda=basisPLDA_org'*dirPLDA2_org;
dirSPCA1_plda=basisPLDA_org'*dirSPCA1_org;
dirSPCA2_plda=basisPLDA_org'*dirSPCA2_org;

Data_spca=basisSPCA_org'*Data_org;
dirORG1_spca=basisSPCA_org'*dirORG1_org;
dirORG2_spca=basisSPCA_org'*dirORG2_org;
dirPCA1_spca=basisSPCA_org'*dirPCA1_org;
dirPCA2_spca=basisSPCA_org'*dirPCA2_org;
dirPLDA1_spca=basisSPCA_org'*dirPLDA1_org;
dirPLDA2_spca=basisSPCA_org'*dirPLDA2_org;
dirSPCA1_spca=basisSPCA_org'*dirSPCA1_org;
dirSPCA2_spca=basisSPCA_org'*dirSPCA2_org;


%%
f1=figure('position',[0 0 900 850],'name',['Linear analysis on data - ' num2str(tag)]); movegui(f1,'south');

subplot 221; hold on; grid on
plot([0 0],[-1 1],'k--'); plot([-1 1],[0 0],'k--');
plot(Data_org(1,find(label==1)),Data_org(2,find(label==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(Data_org(1,find(label==2)),Data_org(2,find(label==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dirORG1_org(1),dirORG1_org(2),'color',CL{1},'linewidth',2,'maxheadsize',.4);
h(2)=quiver(0,0,dirPCA1_org(1),dirPCA1_org(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
h(3)=quiver(0,0,dirPLDA1_org(1),dirPLDA1_org(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
h(4)=quiver(0,0,dirSPCA1_org(1),dirSPCA1_org(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dirORG2_org(1),dirORG2_org(2),'--','color',CL{1},'linewidth',2);
quiver(0,0,dirPCA2_org(1),dirPCA2_org(2),'--','color',CL{2},'linewidth',2);
quiver(0,0,dirPLDA2_org(1),dirPLDA2_org(2),'--','color',CL{3},'linewidth',2);
quiver(0,0,dirSPCA2_org(1),dirSPCA2_org(2),'--','color',CL{4},'linewidth',2);
legend(h(1:4),'original','pca','plda','spca','location','northeast')
xlim([-1 1]); ylim([-1 1]);
title('Original space','interpreter','latex')
set(gca,'fontsize',14); axis square

subplot 222; hold on; grid on
plot([0 0],[-1 1],'k--'); plot([-1 1],[0 0],'k--');
plot(Data_pca(1,find(label==1)),Data_pca(2,find(label==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(Data_pca(1,find(label==2)),Data_pca(2,find(label==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dirORG1_pca(1),dirORG1_pca(2),'color',CL{1},'linewidth',2,'maxheadsize',.4);
h(2)=quiver(0,0,dirPCA1_pca(1),dirPCA1_pca(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
h(3)=quiver(0,0,dirPLDA1_pca(1),dirPLDA1_pca(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
h(4)=quiver(0,0,dirSPCA1_pca(1),dirSPCA1_pca(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dirORG2_pca(1),dirORG2_pca(2),'--','color',CL{1},'linewidth',2);
quiver(0,0,dirPCA2_pca(1),dirPCA2_pca(2),'--','color',CL{2},'linewidth',2);
quiver(0,0,dirPLDA2_pca(1),dirPLDA2_pca(2),'--','color',CL{3},'linewidth',2);
quiver(0,0,dirSPCA2_pca(1),dirSPCA2_pca(2),'--','color',CL{4},'linewidth',2);
% legend(h(1:3),'original','pca','plda','location','northwest')
xlim([-1 1]); ylim([-1 1]);
title('Transformed to PCA space','interpreter','latex')
set(gca,'fontsize',14); axis square

subplot 223; hold on; grid on
plot([0 0],[-1 1],'k--'); plot([-1 1],[0 0],'k--');
plot(Data_plda(1,find(label==1)),Data_plda(2,find(label==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(Data_plda(1,find(label==2)),Data_plda(2,find(label==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dirORG1_plda(1),dirORG1_plda(2),'color',CL{1},'linewidth',2,'maxheadsize',.4);
h(2)=quiver(0,0,dirPCA1_plda(1),dirPCA1_plda(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
h(3)=quiver(0,0,dirPLDA1_plda(1),dirPLDA1_plda(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
h(4)=quiver(0,0,dirSPCA1_plda(1),dirSPCA1_plda(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dirORG2_plda(1),dirORG2_plda(2),'--','color',CL{1},'linewidth',2);
quiver(0,0,dirPCA2_plda(1),dirPCA2_plda(2),'--','color',CL{2},'linewidth',2);
quiver(0,0,dirPLDA2_plda(1),dirPLDA2_plda(2),'--','color',CL{3},'linewidth',2);
quiver(0,0,dirSPCA2_plda(1),dirSPCA2_plda(2),'--','color',CL{4},'linewidth',2);
% legend(h(1:3),'original','pca','plda','location','northwest')
xlim([-1 1]); ylim([-1 1]);
xlabel(['$\alpha$ = ' num2str(alpha)],'interpreter','latex')
title('Transformed to PLDA space','interpreter','latex')
set(gca,'fontsize',14); axis square

subplot 224; hold on; grid on
plot([0 0],[-1 1],'k--'); plot([-1 1],[0 0],'k--');
plot(Data_spca(1,find(label==1)),Data_spca(2,find(label==1)),'o','markerfacecolor',CL{7},'markeredgecolor',CL{7})
plot(Data_spca(1,find(label==2)),Data_spca(2,find(label==2)),'v','markerfacecolor',CL{6},'markeredgecolor',CL{6})
h(1)=quiver(0,0,dirORG1_spca(1),dirORG1_spca(2),'color',CL{1},'linewidth',2,'maxheadsize',.4);
h(2)=quiver(0,0,dirPCA1_spca(1),dirPCA1_spca(2),'color',CL{2},'linewidth',2,'maxheadsize',.4);
h(3)=quiver(0,0,dirPLDA1_spca(1),dirPLDA1_spca(2),'color',CL{3},'linewidth',2,'maxheadsize',.4);
h(4)=quiver(0,0,dirSPCA1_spca(1),dirSPCA1_spca(2),'color',CL{4},'linewidth',2,'maxheadsize',.4);
quiver(0,0,dirORG2_spca(1),dirORG2_spca(2),'--','color',CL{1},'linewidth',2);
quiver(0,0,dirPCA2_spca(1),dirPCA2_spca(2),'--','color',CL{2},'linewidth',2);
quiver(0,0,dirPLDA2_spca(1),dirPLDA2_spca(2),'--','color',CL{3},'linewidth',2);
quiver(0,0,dirSPCA2_spca(1),dirSPCA2_spca(2),'--','color',CL{4},'linewidth',2);
% legend(h(1:3),'original','pca','plda','location','northwest')
xlim([-1 1]); ylim([-1 1]);
xlabel(['$sPCA$ flag = ' num2str(do_spca)],'interpreter','latex')
title(['Transformed to SPCA (Supervised-PCA) space'],'interpreter','latex')
set(gca,'fontsize',14); axis square