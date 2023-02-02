function [cPLDA_directions,cPLDA_projectionsB] = cPLDA_ModesL_ex1v2_cmp_nfe(uB,labelB,ret_dir_num,xxT,lT)
load(['../DATA/METADATA/nfe_params_inside']); % SD_spread=4;
% M=300; N=300;
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
cPLDA_directions=[];
switch alt
    case 1
        AllData=xxT; lN=lT;
        AllData_m0=AllData-repmat(mean(AllData')',1,size(AllData,2));


        [Uall,~,~]=svd(AllData_m0,'econ');
        clear ssall; temp=Uall'*AllData_m0;
        for aa=1:size(temp,1)
            ssall(aa,aa)=1/norm(temp(aa,:));
        end
        for a=1:length(PsiB_m0)
            Psi=PsiB_m0{a};
            PsiPCA=Uall'*Psi;
            PsiPCAB{a}=PsiPCA;
        end

        [Vec,lam] = cPLDA(PsiPCAB,labelB,alpha,ssall); Vec=real(Vec);
        cPLDA_dir=Uall*Vec;


        % %
        [l] = exv2_supp1(cPLDA_dir,PsiB_m0,labelB,SD_spread,lam);
        %         l=1;
        if isempty(l)
            l=1;
        end
        % %


        cPLDA_directions(:,1)=cPLDA_dir(:,l);
        sel_dir=cPLDA_dir(:,l);



        %%


        [B_seldir]=gramschmidt_p1(sel_dir);

        adm0=AllData_m0;

        for b=2:ret_dir_num
            adm0=B_seldir(:,2:end)'*adm0;
            [Uall,~,~]=svd(adm0,'econ');
            clear ssall; temp=Uall'*adm0;
            adm0=B_seldir(:,2:end)*adm0;
            for aa=1:size(temp,1)
                ssall(aa,aa)=1/norm(temp(aa,:));
            end
            for a=1:length(PsiB_m0)
                Psi=PsiB_m0{a};
                Psi=B_seldir(:,2:end)'*Psi;
                PsiPCA=Uall'*Psi;
                PsiPCAB{a}=PsiPCA;
            end



            [Vec,lam] = cPLDA(PsiPCAB,labelB,alpha,ssall); Vec=real(Vec);
            tmp=Uall*Vec;
            cPLDA_dir=B_seldir(:,2:end)*tmp;


            %
            [l] = exv2_supp1(cPLDA_dir,PsiB_m0,labelB,SD_spread,lam);
            %             l=1;
            %

            if isempty(l)
                ret_dir_num=b-1;
                break;
            end
            cPLDA_directions(:,b)=cPLDA_dir(:,l);
            sel_dir=cPLDA_dir(:,l);



            [B_seldir]=gramschmidt_p1(sel_dir);

        end

    otherwise
        disp('ERROR!!');
end

%%
dir_num=ret_dir_num;% min(size(cPLDA_dir,2),2);
% % % cPLDA_directions=cPLDA_dir(:,1:ret_dir_num);
Psi=AllData_m0; Psimean=mean(AllData')';
cPLDA_proj=Psi'*cPLDA_directions;

for a=1:length(PsiB_m0)
    temp=PsiB_m0{a};
    cPLDA_projectionsB{a}=cPLDA_directions'*temp;
end

%%
clear XWC VWC
for a=1:length(labelB)
    label=labelB{a}; cPLDA_projections=cPLDA_projectionsB{a};
    class=sort(unique(label)); pt=15;
    for which=1:size(cPLDA_directions,2)
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

for which=1:size(cPLDA_directions,2)
    if dm(which)>0
    else
        temp=cPLDA_directions(:,which);
        cPLDA_directions(:,which)=-temp;
    end
end

%%
for a=1:length(PsiB_m0)
    temp=PsiB_m0{a};
    cPLDA_projectionsB{a}=cPLDA_directions'*temp;
end
