function [sPLDA_directionsB,sPLDA_projectionsB,viz_pldaB] = sPLDA_ModesL_ex1(uB,labelB,pt_wtB,ret_dir_num,xxT,pT,lT)
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
alt=1; % 1 or 2
does_calc_alpha=0; % 0 means fixed alpha, 1 means calculate alpha by optimization

%%
if does_calc_alpha==0
    alpha=1.618;%8.4237;%4;%size(Psi,2);
elseif does_calc_alpha==1
    % this part under construction
    alpha=[];
    if isempty(alpha)
        Curveoption.low=.01;
        Curveoption.high=10;
        Curveoption.step=0.01;
        Curveoption.nPLDA=3;
        alpha=Calculate_Alpha(Psi,label,EVc,Curveoption);
    end
    alpha
else
    disp('ERROR!!')
end

%%
sPLDA_directionsB=[];
switch alt
    case 1
        %         AllData=xxT; p_wtN=pT; lN=lT;
        %         AllData_m0=AllData-repmat(mean(AllData')',1,size(AllData,2));
        %         [Uall,~,~]=svd(AllData_m0,'econ');

        for a=1:length(PsiB_m0) % % %
            [Uall,~,~]=svd(PsiB_m0{a});

            clear ssall; temp=Uall'*PsiB_m0{a};
            for aa=1:size(temp,1)
                ssall(aa,aa)=1/norm(temp(aa,:));
            end
            %         for a=1:length(PsiB_m0) % % %
            sPLDA_directions=[]; % % %
            Psi=PsiB_m0{a};
            PsiPCA=Uall'*Psi;
            PsiPCAB{a}=PsiPCA;

            [Vec,lam] = PLDA(PsiPCAB{a},labelB{a},alpha); Vec=real(Vec);
            sPLDA_dir=Uall*Vec;

            % % %
            sPLDA_directions(:,1)=sPLDA_dir(:,1);
            sel_dir=sPLDA_dir(:,1);

            [B_seldir]=gramschmidt_p1(sel_dir);

            adm0=PsiB_m0{a};

            for b=2:ret_dir_num
                adm0=B_seldir(:,2:end)'*adm0;
                [Uall,~,~]=svd(adm0); % ,'econ');
                clear ssall; temp=Uall'*adm0;
                adm0=B_seldir(:,2:end)*adm0;
                for aa=1:size(temp,1)
                    ssall(aa,aa)=1/norm(temp(aa,:));
                end

                Psi=PsiB_m0{a};
                Psi=B_seldir(:,2:end)'*Psi;
                PsiPCA=Uall'*Psi;
                PsiPCAB{a}=PsiPCA;

                [Vec,lam] = PLDA(PsiPCAB{a},labelB{a},alpha); Vec=real(Vec);
                tmp=Uall*Vec;
                sPLDA_dir=B_seldir(:,2:end)*tmp;

                sPLDA_directions(:,b)=sPLDA_dir(:,1);
                sel_dir=sPLDA_dir(:,1);
                [B_seldir]=gramschmidt_p1(sel_dir);

            end

            sPLDA_directionsB{a}=sPLDA_directions; % sPLDA_dir(:,1:ret_dir_num);

            % % %
        end

    otherwise
        disp('ERROR!!');
end

%%
dir_num=ret_dir_num; pt_wt=pt_wtB{1}; % p_wtN;

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
