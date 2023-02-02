function [PCA_directions,PCA_projectionsB,viz_plda] = PCA_ModesL(uB,labelB,pt_wtB,ret_dir_num,xxT,pT,lT)
SD_spread=2; M=300; N=300;
for a=1:length(uB)
    u=uB{a};
    [M2,K]=size(u);
    temp=u;
    PsiB{a}=temp;
    Psimean=mean(temp')';
    temp=temp-repmat(Psimean,1,K);
    PsiB_m0{a}=temp;
end

%%
PCA_directions=[];
AllData=xxT; p_wtN=pT; lN=lT;
AllData_m0=AllData-repmat(mean(AllData')',1,size(AllData,2));


[Uall,~,~]=svd(AllData_m0,'econ');
temp=Uall'*AllData_m0;
% for a=1:length(PsiB_m0)
%     Psi=PsiB_m0{a};
%     PsiPCA=Uall'*Psi;
%     PsiPCAB{a}=PsiPCA;
% end

PCA_dir=Uall; %*Vec;
PCA_directions=PCA_dir(:,1:ret_dir_num);

%%
dir_num=ret_dir_num;% min(size(PCA_dir,2),2);
% % % PCA_directions=PCA_dir(:,1:ret_dir_num);
Psi=AllData_m0; Psimean=mean(AllData')'; pt_wt=p_wtN;
PCA_proj=Psi'*PCA_directions;

for i=1:dir_num
    lambda=linspace(-SD_spread*std(PCA_proj(:,i)),SD_spread*std(PCA_proj(:,i)),7);
    big=[];
    for j=1:7
        utemp=Visualize_LOT(Psimean+lambda(j)*PCA_directions(:,i),pt_wt',M,N,2);
        utemp=utemp((M/5)+1:end-(M/5),(M/5)+1:end-(M/5));
        big=[big,mat2gray(utemp)];
    end
    viz_plda(:,:,i)=big;
end
for a=1:length(PsiB_m0)
    temp=PsiB_m0{a};
    PCA_projectionsB{a}=PCA_directions'*temp;
end

%%
clear XWC VWC
for a=1:length(labelB)
    label=labelB{a}; cPLDA_projections=PCA_projectionsB{a};
    class=sort(unique(label)); pt=15;
    for which=1:size(PCA_directions,2)
        for i=1:length(class)
            std_PLDA=SD_spread*std(cPLDA_projections(which,:));
            lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
            [n] =hist(cPLDA_projections(which,label==class(i)),lambdaPLDA');
            XWC(:,i)=n/sum(n);
            VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        end
        mean_values=sum(VWC);
        dm(which,a)=diff(mean_values);
    end
end

dm=mean(dm')';

for which=1:size(PCA_directions,2)
    if dm(which)>0
    else
        temp=PCA_directions(:,which);
        PCA_directions(:,which)=-temp;
        temp=viz_plda(:,:,which);
        viz_plda(:,:,which)=fliplr(temp);
    end
end

%%
for a=1:length(PsiB_m0)
    temp=PsiB_m0{a};
    PCA_projectionsB{a}=PCA_directions'*temp;
end
