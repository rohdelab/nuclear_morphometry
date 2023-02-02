function [fh_all]=proj_n_plot_PLDAl_testC(u_te,label_te,PLDA_directions,viz_plda,which_axis,labelnames)
load(['../DATA/METADATA/params_inside']); % SD_spread=4;

% % %
u_te(isnan(u_te))=0;
% % %

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
fh_all=[];
for ii=which_axis
    %     fh01=figure('position',[-400 40 800 500]);
    fh_all{ii}=figure('position',[-400 40 800 500]);
    % CL={[0.9 0 0],[0.82 0.82 0]}; % CL={[0.27 0.62 0.25],[0.93 0.45 0.10],[0.79 0.15 0.18],[0.09 0.44 0.65]};
    subplot(5,5,1:15); clear XWC VWC WWC gss
    for i=1:length(class)
        std_PLDA=SD_spread*std(PLDA_proj(ii,:));
        lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
        [n] =hist(PLDA_proj(ii,label_te==class(i)),lambdaPLDA');
        XWC(:,i)=n/sum(n);
        VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        WWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]'.*[1:length(XWC(:,i))]';
    end
    EX=sum(VWC); mean_values=EX; MV{ii}=mean_values;
    EX2=sum(WWC); std_values=sqrt(EX2-EX.^2); std_values=sqrt(std_values); SDV{ii}=std_values;
    
    for a=1:length(class)
        temp=gaussian1d(linspace(1,pt,100),std_values(a),mean_values(a));
        gss(:,a)=temp(:);
    end
    
    %     CL={[1 0 1],[0.5 1 0.5],[0 1 1],[0.5 0.5 0.5],[1 1 0],[1 0.5 0.5],[0.5 0.5 1]};
    %     CL={[0.6 0.6 0.6],[0.3 0.3 0.3],[0 0 0],[0.9 0.9 0.9]};
    CL={[0 0 0],[0.3 0.3 0.3],[0.6 0.6 0.6],[0.9 0.9 0.9]};
    for a=1:length(class)
        if a>7
            CL{a}=[rand(1) rand(1) rand(1)];
        end
        bar(-3,0,'facecolor',CL{a}); hold on;
    end
    
    [s,i]=sort(mean_values,'descend');
    ranking=labelnames(i)';
    disp(['Ranking based on ' num2str(which_axis(ii)) ' - th feature:-']); disp(' ')
    disp(ranking);
    
    
    for a=1:length(class)
        %     plot([mean_values(a) mean_values(a)],[0 max(XWC(:))+.08],'--','color',CL{a},'HandleVisibility','off')
        plot(linspace(1,pt,100),gss(:,i(a)),'color',CL{a},'linewidth',2)
        %         area(linspace(1,pt,100),gss(:,i(a)),'facecolor',CL{a},'FaceAlpha',.3)
        area(linspace(1,pt,100),gss(:,i(a)),'facecolor',CL{a},'FaceAlpha',.4,'linewidth',3,'edgecolor',CL{a}); % for powerpoint slide
        %         text(mean_values(i(a)),max(gss(:,i(a)))+.03,labelnames{i(a)},'fontsize',14,'rotation',15)
    end
    xlim([1 pt+0]); ylim([0 max(gss(:))+0.1]);
    set(gca,'XTickLabel',{})
    grid on
    title({['Projection of test data on ' num2str(which_axis(ii)) ' - th cPLDA direction']},'fontsize',15)
    set(gca,'fontsize',14);


    axis_image=viz_plda(:,:,which_axis(ii)); mvim=7;
s=size(axis_image,2)/mvim; st=s/2; xtc0=st:s:size(axis_image,2); xtc=linspace(-SD_spread,SD_spread,mvim);
xtc1=[];
for a=1:length(xtc)
    if a==ceil(mvim/2)
        xtc1{a}=[num2str(xtc(a))];
    else
        xtc1{a}=[num2str(xtc(a)) '\sigma'];
    end
end

    subplot(5,5,[16:25]);
    imshow(viz_plda(:,:,which_axis(ii))); axis on
xticks(xtc0); xticklabels(xtc1); xlabel('Projection score','interpreter','latex'); % xlabel('(\sigma)');
yticks([]);
set(gca,'fontsize',15)
pause(1);
    
    %     fh_all{ii}=fh01;
    
end

mv=0; sdv=0;
for ii=which_axis
    mv=mv+MV{ii};
    sdv=sdv+SDV{ii};
end
mv=mv/length(which_axis); sdv=sdv/length(which_axis);

[s,i]=sort(mv,'descend');
rank_joint=labelnames(i)';
disp(['Joint ranking:-']); disp(' ')
disp(rank_joint);

gss_joint=[];
for a=1:length(class)
    temp=gaussian1d(linspace(1,pt,100),sdv(a),mv(a));
    gss_joint(:,a)=temp(:);
end

fh_all{max(which_axis)+1}=figure('position',[-400 40 800 500]);
subplot(5,5,1:15); hold on
for a=1:length(class)
    plot(linspace(1,pt,100),gss_joint(:,i(a)),'color',CL{a},'linewidth',2)
    area(linspace(1,pt,100),gss_joint(:,i(a)),'facecolor',CL{a},'FaceAlpha',.4,'linewidth',3,'edgecolor',CL{a}); % for powerpoint slide
    %         text(mean_values(i(a)),max(gss(:,i(a)))+.03,labelnames{i(a)},'fontsize',14,'rotation',15)
end
xlim([1 pt+0]); ylim([0 max(gss_joint(:))+0.1]);
set(gca,'XTickLabel',{})
grid on
title({['Projection of test data on all cPLDA direction']},'fontsize',15)
set(gca,'fontsize',14);
subplot(5,5,[16:25]);
imshow(viz_plda(:,:,which_axis(1)),[]);
set(gca,'fontsize',14)

