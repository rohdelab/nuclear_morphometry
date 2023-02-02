function [PLDA_projections,fh01,fh02,mean_values]=proj_n_plot_PLDAl(u_tr,label_tr,u_te,label_te,PLDA_directions,viz_plda,which_axis)
load(['../DATA/METADATA/params_inside']); % SD_spread=4;
M=300; N=300;
[M2_te,K_te]=size(u_te);
Psi=u_te;
Psimean=mean(u_tr')'; %mean(Psi')';
Psi=Psi-repmat(Psimean,1,K_te);

%%
dir_num=length(which_axis);
PLDA_dir_trunc=PLDA_directions(:,which_axis);

class=sort(unique(label_te));
cmap = [1 0 0;0 0 1;0 1 0];
PLDA_proj=PLDA_dir_trunc'*Psi;
pt=15; wdd=1.8;

%%
start_from=1;
fh01=figure('position',[-400 40 800 500]);
% CL={[0.9 0 0],[0.82 0.82 0]};
CL={[0.27 0.62 0.25],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65]};
% CL={[0.1 0.8 0.1],[0.98 0 0],[0.79 0.15 0.18],[0.09 0.44 0.65]}; % colors for powerpoint slide
subplot(5,5,1:15);
for i=1:length(class)
    std_PLDA=SD_spread*std(PLDA_proj(start_from,:));
    lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
    [n] =hist(PLDA_proj(1,label_te==class(i)),lambdaPLDA');
    XWC(:,i)=n/sum(n);
    VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
end
mean_values=sum(VWC);
labelnames1={['Benign (\mu = ' num2str(mean_values(1)) ')'];['Malignant (\mu = ' num2str(mean_values(2)) ')']};
labelnames2={'Benign';'Malignant'}; labelnames1=labelnames2;
for a=1:length(class)
    bar(-3,0,'facecolor',CL{a}); hold on;
end
legend(labelnames1,'location','northwest');
for a=1:length(class)
    plot([mean_values(a) mean_values(a)],[0 max(XWC(:))+.08],'--','color',CL{a},'HandleVisibility','off')
end

for a=1:length(class)
    %     temp=zeros(size(XWC)); temp(:,a)=XWC(:,a);
    %     bar(temp,'facecolor',CL{a}); hold on % consider using histogram
    temp=zeros(length(XWC),3); temp(:,2*a-1)=XWC(:,a);
    bar(temp,wdd,'facecolor',CL{a},'HandleVisibility','off'); hold on % consider using histogram
end
text(mean_values(1),max(XWC(:))+.07,['\Delta\mu = ' num2str(diff(mean_values))],'fontsize',10)
xlim([0 pt+1]);
set(gca,'XTickLabel',{})
grid on
title({['Projection of data on ' num2str(which_axis(1)) ' - th cPLDA direction']},'fontsize',15)
set(gca,'fontsize',10);
subplot(5,5,[16:25]);
imshow(viz_plda(:,:,which_axis(1)),[]);
set(gca,'fontsize',10)
pause(5);

%%
if dir_num>1
    fh02=figure('position',[0 0 1000 800]);movegui(fh02,'northwest')
    set(fh02,'defaultAxesColorOrder',[[0 0 0]; [0 0 0 ]]);
    set(0,'currentfigure',fh02);
    subplot(9,2,[11 13 15])
    which_plda=start_from;
    for a=1:length(class)
        temp=PLDA_proj(which_plda,:);
        std_PLDA1=SD_spread*std(temp);
        lambdaPLDA=linspace(-std_PLDA1,std_PLDA1,pt);
        [n]=hist(temp(label_te==class(a)),lambdaPLDA');
        XWC1(:,a)=n/sum(n);
    end
    for a=1:length(class)
%         temp=zeros(size(XWC1)); temp(:,a)=XWC1(:,a);
%         bar(temp,'facecolor',CL{a}); hold on % consider using histogram
        temp=zeros(length(XWC1),3); temp(:,2*a-1)=XWC1(:,a);
        bar(temp,wdd,'facecolor',CL{a}); hold on % consider using histogram
    end
    xlim([0 pt+1]);
    set(gca,'XTickLabel',{})
    xlabel([num2str(which_axis(1)) ' - th cPLDA direction']);
    grid on
    set(gca,'fontsize',15);
    pause(.3)
    
    
    subplot(2,9,6:8)
    which_plda=start_from+1;
    for a=1:length(class)
        temp=PLDA_proj(which_plda,:);
        std_PLDA2=SD_spread*std(temp);
        lambdaPLDA=linspace(-std_PLDA2,std_PLDA2,pt);
        [n]=hist(temp(label_te==class(a)),lambdaPLDA');
        XWC2(:,a)=n/sum(n);
    end
    for a=1:length(class)
%         temp=zeros(size(XWC2)); temp(:,a)=XWC2(:,a);
%         barh(temp,'facecolor',CL{a}); hold on % consider using histogram
        temp=zeros(length(XWC2),3); temp(:,2*a-1)=XWC2(:,a);
        barh(temp,wdd,'facecolor',CL{a}); hold on % consider using histogram
    end
    ylim([0 pt+1]); set(gca,'YTickLabel',{})
    yyaxis right; set(gca,'YTickLabel',{})
    ylabel([num2str(which_axis(2)) ' - th cPLDA direction']);
    grid on
    set(gca,'fontsize',15);
    pause(.3)
    
    
    subplot 221
    for a=1:length(class)
        temp=PLDA_proj(1:2,:);
        kk{a}=temp(:,label_te==class(a));
    end
    k=kk{1};
    scatter(k(1,:),k(2,:),45,'o','markerfacecolor',CL{1},'markeredgecolor',CL{1},'markerfacealpha',0.3,'markeredgealpha',0.3);hold on
    k=kk{2};
    scatter(k(1,:),k(2,:),45,'o','markerfacecolor',CL{2},'markeredgecolor',CL{2},'markerfacealpha',0.3,'markeredgealpha',0.3);
%     scatter(k(1,:),k(2,:),'ko','markerfacecolor',CL{2});
    if size(kk,2)==3
        k=kk{3};
        scatter(k(1,:),k(2,:),'ko','markerfacecolor','g');
    end
    xlim([-1.1*std_PLDA1 1.1*std_PLDA1]); ylim([-1.1*std_PLDA2 1.1*std_PLDA2]);
    set(gca,'XTick',linspace(-std_PLDA1,std_PLDA1,5))
    set(gca,'XTickLabel',{'-2\sigma','-\sigma','0','\sigma','2\sigma'})
    set(gca,'YTick',linspace(-std_PLDA2,std_PLDA2,5))
    set(gca,'YTickLabel',{'-2\sigma','-\sigma','0','\sigma','2\sigma'})
    grid on;
    ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    title({['Principal cPLDA plane']},'fontsize',12)
    legend(labelnames2,'location','best');
    set(gca,'fontsize',15);
    
    subplot 224
    k=kk{1};
    scatter(k(1,:),k(2,:),'o','markerfacecolor',CL{1},'markeredgecolor',[.45 .45 .45]);hold on
    k=kk{2};
    scatter(k(1,:),k(2,:),'*','markeredgecolor',CL{2});
    if size(kk,2)==3
        k=kk{3};
        scatter(k(1,:),k(2,:),'ko','markerfacecolor','g');
    end
    grid on;
    ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    title({['Whole plane']},'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12);
    
    subplot(9,2,9)
    imshow(viz_plda(:,:,which_axis(1)),[]);
    
    subplot(2,9,5)
    imshow(flipud(viz_plda(:,:,which_axis(2))'),[]);
else
    disp('Number of cPLDA directions too small for a 2D plot.')
end

%%
PLDA_projections=(PLDA_directions'*Psi);
