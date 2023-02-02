function [cPLDA_directions,cPLDA_projectionsB] = scPCA_ModesL_ex1v2_both_cmp_nfe(uB1,labelB1,uB_ref,labelB2_ref,ret_dir_num,xxT1,lT1,xxT_ref,lT_ref)

load(['../DATA/METADATA/nfe_params_inside']); % SD_spread=4;

M=300; N=300;
for a=1:length(uB1)
    u=uB1{a};
    [M2,K]=size(u);
    temp=u;
    PsiB1{a}=temp;
    Psimean1=mean(temp')';
    temp=temp-repmat(Psimean1,1,K);
    PsiB_m0_1{a}=temp;
end

for a=1:length(uB_ref)
    u=uB_ref{a};
    [M2,K]=size(u);
    temp=u;
    PsiB_ref{a}=temp;
    Psimean_ref=mean(temp')';
    temp=temp-repmat(Psimean_ref,1,K);
    PsiB_m0_ref{a}=temp;
end

%%
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

AllData1=xxT1; lN1=lT1;
AllData_m0_1=AllData1-repmat(mean(AllData1')',1,size(AllData1,2));

AllData_ref=xxT_ref; lN_ref=lT_ref;
AllData_m0_ref=AllData_ref-repmat(mean(AllData_ref')',1,size(AllData_ref,2));


[Uall,~,~]=svd(AllData_m0_1,'econ');
clear ssall; temp=Uall'*AllData_m0_1;
for aa=1:size(temp,1)
    ssall(aa,aa)=1/norm(temp(aa,:));
end

% % %

% for a=1:length(PsiB_m0_1)
%     Psi=PsiB_m0_1{a};
%     PsiPCA=Uall'*Psi;
%     PsiPCAB{a}=PsiPCA;
% end

for a=1:length(PsiB_m0_1)
    Psi{a}=PsiB_m0_1{a};
    PsiPCAB{a}=Uall'*Psi{a};
end

% % %

for a=1:length(PsiB_m0_1)
    [Vec,~,~]=svd(PsiPCAB{a}); Vec=real(Vec);

    cPLDA_dir=Uall*Vec;
    %     [l] = exv2_supp1(cPLDA_dir,PsiB_m0_1,labelB1,SD_spread,lam);
    cPLDA_directions{a}(:,1)=cPLDA_dir(:,1); % cPLDA_dir(:,l);
    sel_dir=cPLDA_dir(:,1); % cPLDA_dir(:,l);
    [B_seldir_tmp]=gramschmidt_p1(sel_dir);
    B_seldir{a}=B_seldir_tmp;
    adm0{a}=AllData_m0_1;
end

for ii=1:length(PsiB_m0_1)
    ret_dir_numB{ii}=ret_dir_num;
end

for ii=1:length(PsiB_m0_1)

    for b=2:ret_dir_numB{ii}
        adm0{ii}=B_seldir{ii}(:,2:end)'*adm0{ii};
        [Uall,~,~]=svd(adm0{ii},'econ');
        clear ssall; temp=Uall'*adm0{ii};
        adm0{ii}=B_seldir{ii}(:,2:end)*adm0{ii};
        for aa=1:size(temp,1)
            ssall(aa,aa)=1/norm(temp(aa,:));
        end

        % % %

        %         Psi=PsiB_m0_1{ii};
        %         Psi=B_seldir{ii}(:,2:end)'*Psi;
        %         PsiPCA=Uall'*Psi;

        Psi{ii}=B_seldir{ii}(:,2:end)'*Psi{ii};
        PsiPCA=Uall'*Psi{ii};
        Psi{ii}=B_seldir{ii}(:,2:end)*Psi{ii};

        % % %

        [Vec,~,~]=svd(PsiPCA); Vec=real(Vec);

        tmp=Uall*Vec;
        cPLDA_dir=B_seldir{ii}(:,2:end)*tmp;

        %         [l] = exv2_supp1(cPLDA_dir,PsiB_m0_1,labelB1,SD_spread,lam);
        %         if isempty(l)
        %             ret_dir_numB{ii}=b-1;
        %             break;
        %         end
        cPLDA_directions{ii}(:,b)=cPLDA_dir(:,1); % cPLDA_dir(:,l);
        sel_dir=cPLDA_dir(:,1); % cPLDA_dir(:,l);

        [B_seldir_tmp]=gramschmidt_p1(sel_dir);
        B_seldir{ii}=B_seldir_tmp;
    end

end


%%
Psi=AllData_m0_ref; Psimean=mean(AllData_ref')';

for a=1:length(PsiB_m0_ref)
    temp=PsiB_m0_ref{a};
    cPLDA_projectionsB{a}=cPLDA_directions{a}'*temp;
end

%%
clear XWC VWC
for a=1:length(labelB1)
    label=labelB1{a}; cPLDA_projections=cPLDA_projectionsB{a};
    class=sort(unique(label)); pt=15;
    for which=1:size(cPLDA_directions{a},2)
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

dmB=dm; % mean(dm')';

for a=1:length(labelB1)
    dm=dmB(:,a);
    for which=1:size(cPLDA_directions{a},2)
        if dm(which)>0
        else
            temp=cPLDA_directions{a}(:,which);
            cPLDA_directions{a}(:,which)=-temp;
        end
    end
end

%%
for a=1:length(PsiB_m0_ref)
    temp=PsiB_m0_ref{a};
    cPLDA_projectionsB{a}=cPLDA_directions{a}'*temp;
end
