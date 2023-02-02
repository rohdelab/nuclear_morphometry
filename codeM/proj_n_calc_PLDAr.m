function [PLDA_projections,axis_image,dm]=proj_n_calc_PLDAr(u,label,PLDA_directions,viz_plda)
SD_spread=2;
[M,N,K]=size(u);
Psi=zeros(M*N,K);
for k=1:K
    Psi(:,k)=vec(u(:,:,k));
end
Psimean=mean(Psi')';
Psi=Psi-repmat(Psimean,1,K);

%%
PLDA_projections=PLDA_directions'*Psi;
class=sort(unique(label));
PLDA_proj=PLDA_projections;
pt=15;

%%
axis_image=[];
for which=1:size(PLDA_directions,2)
    for i=1:length(class)
        std_PLDA=SD_spread*std(PLDA_proj(which,:));
        lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
        [n] =hist(PLDA_proj(which,label==class(i)),lambdaPLDA');
        XWC(:,i)=n/sum(n);
        VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
    end
    mean_values=sum(VWC);
    dm(which)=diff(mean_values);
    temp=viz_plda(:,:,which);
    axis_image=[axis_image;temp];
end
