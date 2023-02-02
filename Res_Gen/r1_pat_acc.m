clc
clear all
close all

% delete existing files first inside './CSV_files_pat_acc/' (or later incorporate this part into this code)

%%
tag=1;
subtags={'a','b','c','d'};

% tag=2;
% subtags={'d'};

% tag=3;
% subtags={'a'};

% tag=4;
% subtags={'a'};

%

feat_types={'1','2','3'};

% %

Exp='MAIN'; % 'EXP06_ex1v2'; % 'MAIN'; % 'EXP06_ex1'; 'MAIN_CMP_NFE'  ______ usual - 'MAIN'

single_or_composite=1; % 1 = single, 0 = composite ______ usual - 0 = composite 

% %

clsfrs={'LDA','PLDA','RF','LR','SVM','SVM','kNN'};
kerneltypes={'none','none','none','none','linear','rbf','none'};

clsnmtbl={'LDA';'PLDA';'RF';'LR';'SVM-l';'SVM-k';'kNN'};

%%
p0=pwd;
cd ..
p1=pwd;
cd(p0);

if single_or_composite==1
    res_dir=[p1 '/RESULTS/' Exp '/Dataset_' num2str(tag) '/test01A_PLDA/']
    fnm_str='sPLDA';
else
    res_dir=[p1 '/RESULTS/' Exp '/Dataset_' num2str(tag) '/testA/']
    fnm_str='cPLDA';
end

%%
ur=[]; sr=[]; ue=[]; se=[];
for a=1:length(clsfrs)
    u_tr_B=[]; s_tr_B=[]; u_te_B=[]; s_te_B=[];
    for b=1:length(feat_types)
        nm=['pat_acc_feat' feat_types{b} '_' clsfrs{a} kerneltypes{a} '.mat'];
        load([res_dir nm]);
        u_tr_B=[u_tr_B u_tr];
        s_tr_B=[s_tr_B s_tr];
        u_te_B=[u_te_B u_te];
        s_te_B=[s_te_B s_te];
    end
    ur=[ur;u_tr_B];
    sr=[sr;s_tr_B];
    ue=[ue;u_te_B];
    se=[se;s_te_B];
end


%%
clear tr te
for a=1:size(ur,1)
    for b=1:size(ur,2)
        tr{a,b}=[num2str(ur(a,b),'%.0f') ' $\pm$ ' num2str(sr(a,b),'%.1f')];
        te{a,b}=[num2str(ue(a,b),'%.0f') ' $\pm$ ' num2str(se(a,b),'%.1f')];
    end
end

tr
te


%%
% load patients.mat
% T = table(LastName,Age,Weight,Smoker);
% T(1:5,:)

% T_tr=table(tr); T_te=table(te);
% filename_tr = './CSV_files/pat_acc_tr.xlsx';
% filename_te = './CSV_files/pat_acc_te.xlsx';
% writetable(T_tr,filename_tr,'Sheet',1);
% writetable(T_te,filename_te,'Sheet',1);



Liver1=tr(:,1); Thyroid1=tr(:,2); Mesothelioma1=tr(:,3); Melanoma1=tr(:,4);
Liver2=tr(:,5); Thyroid2=tr(:,6); Mesothelioma2=tr(:,7); Melanoma2=tr(:,8);
Liver3=tr(:,9); Thyroid3=tr(:,10); Mesothelioma3=tr(:,11); Melanoma3=tr(:,12);
T_tr=table(clsnmtbl,Liver1,Thyroid1,Mesothelioma1,Melanoma1,...
    Liver2,Thyroid2,Mesothelioma2,Melanoma2,...
    Liver3,Thyroid3,Mesothelioma3,Melanoma3);
Liver1=te(:,1); Thyroid1=te(:,2); Mesothelioma1=te(:,3); Melanoma1=te(:,4);
Liver2=te(:,5); Thyroid2=te(:,6); Mesothelioma2=te(:,7); Melanoma2=te(:,8);
Liver3=te(:,9); Thyroid3=te(:,10); Mesothelioma3=te(:,11); Melanoma3=te(:,12);
T_te=table(clsnmtbl,Liver1,Thyroid1,Mesothelioma1,Melanoma1,...
    Liver2,Thyroid2,Mesothelioma2,Melanoma2,...
    Liver3,Thyroid3,Mesothelioma3,Melanoma3);
filename_tr = ['./CSV_files_pat_acc/' fnm_str '_pat_acc_tr.xlsx'];
filename_te = ['./CSV_files_pat_acc/' fnm_str '_pat_acc_te.xlsx'];
writetable(T_tr,filename_tr,'Sheet',1);
writetable(T_te,filename_te,'Sheet',1);




