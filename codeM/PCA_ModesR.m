function [PCA_directions,PCA_projections,viz_pca] = PCA_ModesR(u,I0,ret_dir_num)
SD=1;
[M,N,K]=size(u);
Psi=zeros(M*N,K);
for k=1:K
    Psi(:,k)=vec(u(:,:,k));
end
Psimean=mean(Psi')';
Psi=Psi-repmat(Psimean,1,K);
[U,S,V]=svd(Psi,'econ');

PCA_dir=U;
dir_num=3;

viz_pca=[];
viz_pca_rcdt=[];
stdev=1.5*mean(std(PCA_dir'*Psi));
for i=1:dir_num
    big=[];
    bigRCDT=[];
    lambda=linspace(-SD*max(std(PCA_dir(:,i)'*Psi),stdev),SD*max(std(PCA_dir(:,i)'*Psi),stdev),5);
    for j=1:5
        utemp=reshape(Psimean+lambda(j)*PCA_dir(:,i),M,N);
        bigRCDT=[bigRCDT,utemp];
        big=[big,iRCDT(utemp,I0)];
    end
    viz_pca=[viz_pca;big];
    viz_pca_rcdt=[viz_pca_rcdt;bigRCDT];
end
PCA_directions=PCA_dir(:,1:ret_dir_num);
PCA_projections=(PCA_directions'*Psi);

