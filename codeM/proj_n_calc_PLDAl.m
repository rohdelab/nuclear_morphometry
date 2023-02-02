function [PLDA_projections,axis_image,dm,dt,dt_tr]=proj_n_calc_PLDAl(u,label,PLDA_directions,viz_plda,projPLDA_tr,lab_tr)
load(['../DATA/METADATA/params_inside']); % SD_spread=4;
M=300; N=300;
[M2,K]=size(u);
Psi=u;
Psimean=mean(Psi')';
Psi=Psi-repmat(Psimean,1,K);

%%
PLDA_projections=PLDA_directions'*Psi;
class=sort(unique(label));
PLDA_proj=PLDA_projections;
pt=15;
pdf_pt=100;
eps=1e-8; rm_edge = 1;



%%
axis_image=[];
for which=1:size(PLDA_directions,2)
    for i=1:length(class)
        std_PLDA=SD_spread*std(PLDA_proj(which,:));
        lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
        [n] =hist(PLDA_proj(which,label==class(i)),lambdaPLDA');
        XWC(:,i)=n/sum(n);
        VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        
        tmp=spline(1:pt,XWC(:,i),linspace(1,pt,pdf_pt));
        pdf_XWC(:,i)=tmp/sum(tmp);

        I = pdf_XWC(:,i);
        I = abs(I)+eps;
        I_domain = linspace(0,1,length(I));    % domain of 1D signal
        Ihat_domain = linspace(0,1,length(I)); % domain of CDT
        Ihat = CDT2(I_domain, I/sum(I), Ihat_domain, rm_edge); % CDT of each sample
        
        cdt_XWC(:,which,i)=Ihat;
    end
    
%     %
%     f=figure('position',[0 0 500 600]); movegui(f,'north');
%     subplot 311
%     bar(XWC); grid on; set(gca,'fontsize',12)
%     subplot 312
%     plot(pdf_XWC); grid on; set(gca,'fontsize',12)
%     subplot 313
%     plot(squeeze(cdt_XWC(:,which,:))); grid on; set(gca,'fontsize',12)
%     %
    
    for i=1:length(class)
        std_PLDA=SD_spread*std(projPLDA_tr(which,:));
        lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
        [n] =hist(projPLDA_tr(which,lab_tr==class(i)),lambdaPLDA');
        XWC_tr(:,i)=n/sum(n);
        VWC_tr(:,i)=XWC_tr(:,i).*[1:length(XWC_tr(:,i))]';
        
        tmp=spline(1:pt,XWC_tr(:,i),linspace(1,pt,pdf_pt));
        pdf_XWC_tr(:,i)=tmp/sum(tmp);

        I = pdf_XWC_tr(:,i);
        I = abs(I)+eps;
        I_domain = linspace(0,1,length(I));    % domain of 1D signal
        Ihat_domain = linspace(0,1,length(I)); % domain of CDT
        Ihat = CDT2(I_domain, I/sum(I), Ihat_domain, rm_edge); % CDT of each sample
        
        cdt_XWC_tr(:,which,i)=Ihat;
    end
    
%     %
%     f=figure('position',[0 0 500 600]); movegui(f,'north');
%     subplot 311
%     bar(XWC_tr); grid on; set(gca,'fontsize',12)
%     subplot 312
%     plot(pdf_XWC_tr); grid on; set(gca,'fontsize',12)
%     subplot 313
%     plot(squeeze(cdt_XWC_tr(:,which,:))); grid on; set(gca,'fontsize',12)
%     %
    
    
    
    
    mean_values=sum(VWC);
    dm(which)=diff(mean_values);
    temp=viz_plda(:,:,which);
    axis_image=[axis_image;temp];
end

for i=1:length(class)
    temp=cdt_XWC(:,:,i);
    temp_tr=cdt_XWC_tr(:,:,i);
    s(:,i)=temp(:);
    s_tr(:,i)=temp_tr(:);
end
dt=norm(s(:,1)-s(:,2));
dt_tr=(norm(s(:,1)-s_tr(:,1))+norm(s(:,2)-s_tr(:,2)))/2;
