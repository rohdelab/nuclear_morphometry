function PCA_projections = proj_n_plot_PCAr(u,PCA_directions,viz_pca)
SD=1;
[M,N,K]=size(u);
Psi=zeros(M*N,K);
for k=1:K
    Psi(:,k)=vec(u(:,:,k));
end
Psimean=mean(Psi')';
Psi=Psi-repmat(Psimean,1,K);

%%
dir_num=2;
PCA_dir_trunc=PCA_directions(:,1:dir_num);

PCA_proj=(PCA_dir_trunc'*Psi)';

%%
fh1=figure('position',[0 0 650 300]);
subplot(121)
imshow(viz_pca,[])
title('Inverse of R-CDT PCA modes')
xlabel('Variation along a mode');
ylabel('PCA modes');
set(gca,'fontsize',14);

subplot(1,2,2)
scatter(PCA_proj(:,1),PCA_proj(:,2),24,'b','markerfacecolor','k')
grid on; title('Principal PCA plane');
ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
set(gca,'fontsize',14);



% for a=1:length(class)
%     temp=PCA_proj(1:2,:);
%     kk{a}=temp(:,label==class(a));
% end
% k=kk{1};
% scatter(k(1,:),k(2,:),[CL{1} '*']);hold on
% k=kk{2};
% scatter(k(1,:),k(2,:),'ko','markerfacecolor',CL{2});
% grid on;
% ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
% title({['PLDA plane']},'fontsize',12)
% legend({'Benign','Malignant'},'location','best');
% set(gca,'fontsize',15);



PCA_projections=(PCA_directions'*Psi);
