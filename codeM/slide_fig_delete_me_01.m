clc
clear all
close all

warning('off','all')

%%
% tag=0; subtag=[]; nme='toy'; ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'}; % toy data % subtag is a character like 'a', 'b', etc.

% tag=1; subtag='a'; nme='liver'; ote={'1b','1c','1d'}; otnme={'thyroid','mesothelioma','melanoma'}; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; nme='thyroid'; ote={'1a','1c','1d'}; otnme={'liver','mesothelioma','melanoma'}; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; nme='mesothelioma'; ote={'1a','1b','1d'}; otnme={'liver','thyroid','melanoma'}; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
tag=1; subtag='d'; nme='melanoma'; ote={'1a','1b','1c'}; otnme={'liver','thyroid','mesothelioma'}; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

%%
p0=pwd; cd ..; pp=pwd;
mdpth=[pwd '/DATA/METADATA'];
respth=[pwd '/RESULTS/EXP01'];

cd(p0);

%%
ret_dir_num=20; Nfold=10;

foldset=5;
vnm1=['Exp01_data' num2str(tag) subtag '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(ret_dir_num) '_PLDA_directions'];
vnm2=['Exp01_data' num2str(tag) subtag '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(ret_dir_num) '_projPLDA_tr'];
vnm3=['Exp01_data' num2str(tag) subtag '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(ret_dir_num) '_viz_plda'];


vnm1pnm=[mdpth '/EXP01/' vnm1];
vnm2pnm=[mdpth '/EXP01/' vnm2];
vnm3pnm=[mdpth '/EXP01/' vnm3];



disp('Calculating and saving...');
fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 10-fold split

vnm4=['Exp01_data' num2str(tag) subtag num2str(tag) subtag '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(ret_dir_num) '_projPLDA_te'];
vnm4pnm=[mdpth '/EXP01/' vnm4];

indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
indpnm=[mdpth '/' indnm];

load(indpnm)
load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);

inp_te=[pp '/DATA/data' num2str(tag) subtag '/rcdt'];
load([inp_te '/Rcdt']);
data_te=u; label_te=label;

[~,~,xx_te,lab_te]=fn_trte_split(data_te,label_te,ind,fold_tr,fold_te);
[projPLDA_te,f01,f02]=proj_n_plot_PLDAr(xx_te,lab_te,PLDA_directions,viz_plda);

% 
% for a=1:length(ote)
%     vnm4=['Exp01_data' num2str(tag) subtag ote{a} '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(ret_dir_num) '_projPLDA_te'];
%     vnm4pnm=[mdpth '/EXP01/' vnm4];
%     
%     indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
%     indpnm=[mdpth '/' indnm];
%     
%     load(indpnm)
%     load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
%     
%     inp_te=[pp '/DATA/data' ote{a} '/rcdt'];
%     load([inp_te '/Rcdt']);
%     data_te=u; label_te=label;
%     xx_te=data_te;  lab_te=label_te;
%     [projPLDA_te,f01,f02]=proj_n_plot_PLDAr(xx_te,lab_te,PLDA_directions,viz_plda);
% end


