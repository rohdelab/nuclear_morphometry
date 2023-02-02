function [PLDA_projections,fh01,fh02,mean_values]=proj_n_plot_PLDAl_ResGen(u_tr,label_tr,u_te,label_te,PLDA_directions,viz_plda,which_axis,szmltplr,dtnm)
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
pt=10; wdd=1.8;

%%
start_from=1;
fh01=figure('position',[-400 40 400*szmltplr 250*szmltplr]); movegui(fh01,'northwest')
% CL={[0.9 0 0],[0.82 0.82 0]};
% CL={[0.27 0.62 0.25],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65]};
CL={[0.27 0.62 0.25],[0.9 0.7 0.12],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65]};
CL={[0.27 0.62 0.25],[0.79 0.15 0.18],[0.09 0.44 0.65]};
% CL={[0.1 0.8 0.1],[0.98 0 0],[0.79 0.15 0.18],[0.09 0.44 0.65]}; % colors for powerpoint slide
subplot(5,5,6:20);
for i=1:length(class)
    std_PLDA=SD_spread*std(PLDA_proj(start_from,:));
    lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
    [n] =hist(PLDA_proj(1,label_te==class(i)),lambdaPLDA');
    XWC(:,i)=n/sum(n);
    VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
end
mean_values=sum(VWC);
labelnames1={['Benign (\mu = ' num2str(mean_values(1)) ')'];['Malignant (\mu = ' num2str(mean_values(2)) ')']};
labelnames2={'Less (or not) malig';'More malig'}; labelnames1=labelnames2;
for a=1:length(class)
    bar(-3,0,'facecolor',CL{a},'edgecolor',CL{a}); hold on;
end
legend(labelnames1,'location','northeast'); %,'interpreter','latex');
legend boxoff
for a=1:length(class)
    plot([mean_values(a) mean_values(a)],[0 max(XWC(:))+.08],'--','linewidth',2,'color',CL{a},'HandleVisibility','off')
end

for a=1:length(class)
    %     temp=zeros(size(XWC)); temp(:,a)=XWC(:,a);
    %     bar(temp,'facecolor',CL{a}); hold on % consider using histogram
    temp=zeros(length(XWC),3); temp(:,2*a-1)=XWC(:,a);
    bar(temp,wdd,'facecolor',CL{a},'edgecolor',CL{a},'HandleVisibility','off'); hold on % consider using histogram
end
% text(mean_values(1),max(XWC(:))+.07,['\Delta\mu = ' num2str(diff(mean_values))],'fontsize',10)
% text(mean_values(1)+diff(mean_values)/4,max(XWC(:))+.07,['\Delta\mu'],'fontsize',15)
text(mean_values(1)-1,max(XWC(:))+.05,['\Delta\mu'],'fontsize',15)
ylabel('Percentage incidents'); %,'interpreter','latex');
% ylabel('Projection of test data'); %,'interpreter','latex');
xlim([0 pt+1]);
set(gca,'XTickLabel',{})
% grid on
axis off
% title({['Projection of ' dtnm ' test data']},'fontsize',18,'interpreter','latex')
set(gca,'fontsize',15);
pause(1)

axis_image=viz_plda(:,:,which_axis(1)); mvim=7;
s=size(axis_image,2)/mvim; st=s/2;
xtc0=st:s:size(axis_image,2); xtc=linspace(-SD_spread,SD_spread,mvim);

xtc1=[];
for a=1:length(xtc)
    if a==ceil(mvim/2)
        xtc1{a}=[num2str(xtc(a))];
    else
        xtc1{a}=[num2str(xtc(a)) '\sigma'];
    end
end

fh02=figure('position',[-400 40 405*szmltplr 180*szmltplr]); movegui(fh02,'northeast')
subplot(2,5,[1:5]);

axis_image=viz_plda(:,:,which_axis(1));
% %
% [axis_image] = ib2w(axis_image,0);
% %
imshow(axis_image); axis on;

% imagesc(viz_plda(:,:,which_axis(1))); colormap gray
xticks(xtc0); xticklabels(xtc1); xlabel('Projection score','interpreter','latex'); % xlabel('(\sigma)');
yticks([]);
set(gca,'fontsize',15)
pause(1);

%%
PLDA_projections=(PLDA_directions'*Psi);
