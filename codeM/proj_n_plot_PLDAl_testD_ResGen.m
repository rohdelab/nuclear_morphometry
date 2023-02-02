function [fh_all1,fh_all2]=proj_n_plot_PLDAl_testD_ResGen(u_te,label_te,PLDA_directions,viz_plda,which_axis,labelnames,szmltplr,dtnm)
load(['../DATA/METADATA/params_inside']); % SD_spread=4;
M=300; N=300;
[M2_te,K_te]=size(u_te);
Psi=u_te;
Psimean=mean(u_te')'; %mean(Psi')';
Psi=Psi-repmat(Psimean,1,K_te);

%%
dir_num=length(which_axis);
PLDA_dir_trunc=PLDA_directions(:,which_axis);

class=sort(unique(label_te));
cmap = [1 0 0;0 0 1;0 1 0];
PLDA_proj=PLDA_dir_trunc'*Psi;
pt=10; wdd=1.8;

%%
% CL={[0.27 0.62 0.25],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65
% CL={[0.27 0.62 0.25],[0.9 0.7 0.12],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65]};
CL={[0.27 0.62 0.25],[0.79 0.15 0.18],[0.09 0.44 0.65]};
fh_all1=[]; fh_all2=[];
for ii=which_axis
    fh_all1{ii}=figure('position',[-400 40 400*szmltplr 250*szmltplr]); movegui(fh_all1{ii},'northwest')
    subplot(5,5,6:20); clear XWC VWC WWC
    for i=1:length(class)
        std_PLDA=SD_spread*std(PLDA_proj(ii,:));
        lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
        [n] =hist(PLDA_proj(ii,label_te==class(i)),lambdaPLDA');
        XWC(:,i)=n/sum(n);
        VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        WWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]'.*[1:length(XWC(:,i))]';
    end
    EX=sum(VWC); mean_values=EX;
    EX2=sum(WWC); std_values=sqrt(EX2-EX.^2); std_values=sqrt(std_values);

    mean_values=sum(VWC);
    for a=1:length(class)
        bar(-3,0,'facecolor',CL{a},'edgecolor',CL{a}); hold on;
    end
    legend(labelnames,'location','northwest'); legend boxoff
    for a=1:length(class)
        plot([mean_values(a) mean_values(a)],[0 max(XWC(:))+.08],'--','linewidth',3,'color',CL{a},'HandleVisibility','off')
    end

    for a=1:length(class)
        temp=zeros(length(XWC),3); temp(:,2*a-1)=XWC(:,a);
        bar(temp,wdd,'facecolor',CL{a},'edgecolor',CL{a},'HandleVisibility','off'); hold on % consider using histogram
    end
    %     text(mean_values(1),max(XWC(:))+.07,['\Delta\mu = ' num2str(diff(mean_values))],'fontsize',10)
    xlim([0 pt+1]);
    set(gca,'XTickLabel',{})
    %     grid on
    axis off
    %     title({['Projection of ' dtnm ' test data']},'fontsize',18,'interpreter','latex')
    set(gca,'fontsize',18);
    pause(1)


    axis_image=viz_plda(:,:,which_axis(1)); mvim=7;
    s=size(axis_image,2)/mvim; st=s/2; xtc0=st:s:size(axis_image,2); xtc=linspace(-SD_spread,SD_spread,mvim);

    xtc1=[];
    for a=1:length(xtc)
        if a==ceil(mvim/2)
            xtc1{a}=[num2str(xtc(a))];
        else
            xtc1{a}=[num2str(xtc(a)) '\sigma'];
        end
    end

    fh_all2{ii}=figure('position',[-400 40 405*szmltplr 180*szmltplr]); movegui(fh_all2{ii},'northeast')
    subplot(2,5,[1:5]);

    %     imagesc(viz_plda(:,:,which_axis(ii))); colormap gray

    axis_image=viz_plda(:,:,which_axis(ii));
    % %
%     [axis_image] = ib2w(axis_image,0);
    % %

    imshow(axis_image); axis on;

    xticks(xtc0); xticklabels(xtc1); xlabel('Projection score','interpreter','latex'); % xlabel('(\sigma)');
    yticks([]);
    set(gca,'fontsize',18)

    [s,i]=sort(mean_values,'descend');
    ranking=labelnames(i)';
    disp(['Ranking based on ' num2str(which_axis(ii)) ' - th feature:-']); disp(' ')
    disp(ranking);


    pause(1);

end
