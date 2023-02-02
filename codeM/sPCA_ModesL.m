function [sPLDA_directionsB,sPLDA_projectionsB,viz_pldaB] = sPCA_ModesL(uB,labelB,pt_wtB,ret_dir_num,xxT,pT,lT)
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
sPLDA_directionsB=[];

AllData=xxT; p_wtN=pT; lN=lT;
AllData_m0=AllData-repmat(mean(AllData')',1,size(AllData,2));

[Uall,~,~]=svd(AllData_m0,'econ');
temp=Uall'*AllData_m0;

for a=1:length(PsiB_m0)
    Psi=PsiB_m0{a};

    [Vec,~,~]=svd(Psi); Vec=real(Vec);
    sPLDA_dir=Vec;

    sPLDA_directionsB{a}=sPLDA_dir(:,1:ret_dir_num);
end

%%
dir_num=ret_dir_num; pt_wt=p_wtN;

viz_pldaB=[];
for a=1:length(PsiB_m0)

    Psi=PsiB_m0{a}; Psimean=mean(PsiB{a}')';

    sPLDA_directions=sPLDA_directionsB{a};
    sPLDA_proj=Psi'*sPLDA_directions;

    viz_plda=[];
    for i=1:dir_num
        lambda=linspace(-SD_spread*std(sPLDA_proj(:,i)),SD_spread*std(sPLDA_proj(:,i)),7);
        big=[];
        for j=1:7
            utemp=Visualize_LOT(Psimean+lambda(j)*sPLDA_directions(:,i),pt_wt',M,N,2);
            utemp=utemp((M/5)+1:end-(M/5),(M/5)+1:end-(M/5));
            big=[big,mat2gray(utemp)];
        end
        viz_plda(:,:,i)=big;
    end

    temp=PsiB_m0{a};
    sPLDA_projectionsB{a}=sPLDA_directions'*temp;
    viz_pldaB{a}=viz_plda;
end
%%
clear XWC VWC
for a=1:length(labelB)
    label=labelB{a}; sPLDA_projections=sPLDA_projectionsB{a};
    class=sort(unique(label)); pt=15;
    for which=1:dir_num %size(sPLDA_directions,2)
        for i=1:length(class)
            std_PLDA=SD_spread*std(sPLDA_projections(which,:));
            lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
            [n] =hist(sPLDA_projections(which,label==class(i)),lambdaPLDA');
            XWC(:,i)=n/sum(n);
            VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        end
        mean_values=sum(VWC);
        dm(which,a)=diff(mean_values);
    end
end

% dm=mean(dm')';
dmB=dm;

for a=1:length(PsiB_m0)
    dm=dmB(:,a);
    sPLDA_directions=sPLDA_directionsB{a};
    viz_plda=viz_pldaB{a};
    for which=1:dir_num % size(sPLDA_directions,2)
        if dm(which)>0
        else
            temp=sPLDA_directions(:,which);
            sPLDA_directions(:,which)=-temp;
            temp=viz_plda(:,:,which);
            viz_plda(:,:,which)=fliplr(temp);
        end
    end

    sPLDA_directionsB{a}=sPLDA_directions;
    viz_pldaB{a}=viz_plda;
    temp=PsiB_m0{a};
    sPLDA_projectionsB{a}=sPLDA_directions'*temp;
end
