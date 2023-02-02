function [PLDA_projections,fh01,fh02,mean_values]=proj_n_plot_PLDAr(u,label,PLDA_directions,viz_plda,which_axis)
SD_spread=2; 
[M,N,K]=size(u);
Psi=zeros(M*N,K);
for k=1:K
    Psi(:,k)=vec(u(:,:,k));
end
Psimean=mean(Psi')';
Psi=Psi-repmat(Psimean,1,K);

%%
dir_num=min(size(PLDA_directions,2),2);
PLDA_dir_trunc=PLDA_directions(:,which_axis:which_axis+dir_num-1);
class=sort(unique(label));
cmap = [1 0 0;0 0 1;0 1 0];
PLDA_proj=PLDA_dir_trunc'*Psi;
pt=15; wdd=1.8;

%%
which=1;
fh01=figure('position',[-400 40 800 500]);
CL={[0.9 0 0],[0.82 0.82 0]};
subplot(5,5,1:15);
for i=1:length(class)
    std_PLDA=SD_spread*std(PLDA_proj(which,:));
    lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
    [n] =hist(PLDA_proj(1,label==class(i)),lambdaPLDA');
    XWC(:,i)=n/sum(n);
    VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
end
mean_values=sum(VWC);
labelnames1={['Benign (\mu = ' num2str(mean_values(1)) ')'];['Malignant (\mu = ' num2str(mean_values(2)) ')']};
labelnames2={'Benign';'Malignant'};
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
title({'Projection of data on the first PLDA direction '},'fontsize',15)
set(gca,'fontsize',10);
subplot(5,5,[16:25]);
imshow(viz_plda(:,:,which_axis),[]);
set(gca,'fontsize',10)
pause(5);

%%
if dir_num>1
    fh02=figure('position',[0 0 1000 800]);movegui(fh02,'northwest')
    set(fh02,'defaultAxesColorOrder',[[0 0 0]; [0 0 0 ]]);
    set(0,'currentfigure',fh02);
    subplot(9,2,[11 13 15])
    which_plda=which;
    for a=1:length(class)
        temp=PLDA_proj(which_plda,:);
        std_PLDA1=SD_spread*std(temp);
        lambdaPLDA=linspace(-std_PLDA1,std_PLDA1,pt);
        [n]=hist(temp(label==class(a)),lambdaPLDA');
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
    xlabel([num2str(which_plda) '^{st} PLDA direction']);
    grid on
    set(gca,'fontsize',15);
    pause(.3)
    
    
    subplot(2,9,6:8)
    which_plda=which+1;
    for a=1:length(class)
        temp=PLDA_proj(which_plda,:);
        std_PLDA2=SD_spread*std(temp);
        lambdaPLDA=linspace(-std_PLDA2,std_PLDA2,pt);
        [n]=hist(temp(label==class(a)),lambdaPLDA');
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
    ylabel([num2str(which_plda) '^{nd} PLDA direction']);
    grid on
    set(gca,'fontsize',15);
    pause(.3)
    
    
    subplot 221
    for a=1:length(class)
        temp=PLDA_proj(1:2,:);
        kk{a}=temp(:,label==class(a));
    end
    k=kk{1};
    scatter(k(1,:),k(2,:),'o','markerfacecolor',CL{1},'markeredgecolor',[.45 .45 .45]);hold on
    k=kk{2};
    scatter(k(1,:),k(2,:),'*','markeredgecolor',CL{2});
%     scatter(k(1,:),k(2,:),'ko','markerfacecolor',CL{2});
    if size(kk,2)==3
        k=kk{3};
        scatter(k(1,:),k(2,:),'ko','markerfacecolor','g');
    end
    xlim([-1.1*std_PLDA1 1.1*std_PLDA1]); ylim([-1.1*std_PLDA2 1.1*std_PLDA2]);
    set(gca,'XTick',linspace(-std_PLDA1,std_PLDA1,5))
    set(gca,'XTickLabel',{'-2\sigma','\sigma','0','\sigma','2\sigma'})
    set(gca,'YTick',linspace(-std_PLDA2,std_PLDA2,5))
    set(gca,'YTickLabel',{'-2\sigma','\sigma','0','\sigma','2\sigma'})
    grid on;
    ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    title({['Principal PLDA plane']},'fontsize',12)
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
    imshow(viz_plda(:,:,which_axis),[]);
    
    subplot(2,9,5)
    imshow(viz_plda(:,:,which_axis+1)',[]);
else
    disp('Number of PLDA directions too small for a 2D plot.')
end

%%
PLDA_projections=(PLDA_directions'*Psi);
