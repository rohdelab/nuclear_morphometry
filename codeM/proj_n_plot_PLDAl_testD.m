function [fh_all]=proj_n_plot_PLDAl_testD(u_te,label_te,PLDA_directions,viz_plda,which_axis,labelnames)
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
pt=15; wdd=1.8;

%%
CL={[0.27 0.62 0.25],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65]};
fh_all=[];
for ii=which_axis
    fh_all{ii}=figure('position',[-400 40 800 500]);
    subplot(5,5,1:15); clear XWC VWC WWC 
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
        bar(-3,0,'facecolor',CL{a}); hold on;
    end
    legend(labelnames,'location','northwest');
    for a=1:length(class)
        plot([mean_values(a) mean_values(a)],[0 max(XWC(:))+.08],'--','linewidth',3,'color',CL{a},'HandleVisibility','off')
    end
    
    for a=1:length(class)
        temp=zeros(length(XWC),3); temp(:,2*a-1)=XWC(:,a);
        bar(temp,wdd,'facecolor',CL{a},'HandleVisibility','off'); hold on % consider using histogram
    end
%     text(mean_values(1),max(XWC(:))+.07,['\Delta\mu = ' num2str(diff(mean_values))],'fontsize',10)
    xlim([0 pt+1]);
    set(gca,'XTickLabel',{})
    grid on
    title({['Projection of data on ' num2str(which_axis(ii)) ' - th cPLDA direction']},'fontsize',15)
    set(gca,'fontsize',22);
    subplot(5,5,[16:25]);
    imshow(viz_plda(:,:,which_axis(ii)),[]);
    
    [s,i]=sort(mean_values,'descend');
    ranking=labelnames(i)';
    disp(['Ranking based on ' num2str(which_axis(ii)) ' - th feature:-']); disp(' ')
    disp(ranking);
    
    
    pause(5);
    
    %     fh_all{ii}=fh01;
    
end
