function [PLDA_projections]=proj_n_plot_PLDAl_cmp_nfe(u_tr,label_tr,u_te,label_te,PLDA_directions,which_axis)
load(['../DATA/METADATA/nfe_params_inside']); % SD_spread=4;
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

for i=1:length(class)
    std_PLDA=SD_spread*std(PLDA_proj(start_from,:));
    lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
    [n] =hist(PLDA_proj(1,label_te==class(i)),lambdaPLDA');
    XWC(:,i)=n/sum(n);
    VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
end
mean_values=sum(VWC);

%%
PLDA_projections=(PLDA_directions'*Psi);
