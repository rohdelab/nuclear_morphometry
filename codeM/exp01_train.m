clc
clear all
close all

warning('off','all')

%% DATA SELECT
tag=0; subtag=[]; % toy data % subtag is a character like 'a', 'b', etc.

% tag=1; subtag='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

%% PARAMETER SELECT
Ndir=20; Nfold=10;
force_run=1;% 1 for must run no matter what, 0 for may skip if already there.
reg_str='TOF';

%% LOAD DATA & PRELIMINARIES 
p0=pwd; cd ..
inp_tr=[pwd '/DATA/data' num2str(tag) subtag '/rcdt'];
mdpth=[pwd '/DATA/METADATA'];
logpth=[pwd '/RESULTS/EXP01/log'];

load([inp_tr '/Rcdt_' reg_str]);
data_tr=u; I0_tr=I0; label_tr=label;
Exp_str = 'Exp01'; Exp_strF = '/EXP01/'; 
cd(p0);

%% CALCULATION
wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
for foldset=1:Nfold
    
    vnm1=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
    vnm2=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
    vnm3=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
    
    vnm1pnm=[mdpth Exp_strF vnm1]; 
    vnm2pnm=[mdpth Exp_strF vnm2]; 
    vnm3pnm=[mdpth Exp_strF vnm3];
    
    ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
    
    if force_run|ISthere1==0|ISthere2==0|ISthere3==0
        disp([Exp_str '-train on Data ' num2str(tag) subtag ': Calculating and saving... ' num2str(foldset) ' of ' num2str(Nfold)]);
        
        fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 10-fold split
        
        indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
        indpnm=[mdpth '/' indnm];
        
        load(indpnm)
        
        [xx_tr,lab_tr,~,~]=fn_trte_split(data_tr,label_tr,ind,fold_tr,fold_te);
        [PLDA_directions,projPLDA_tr,viz_plda]=PLDA_ModesR(xx_tr,lab_tr,I0_tr,1e+4,Ndir);
        
        
        save(vnm1pnm,'PLDA_directions','-v7.3');
        save(vnm2pnm,'projPLDA_tr','-v7.3');
        save(vnm3pnm,'viz_plda','-v7.3');
        pause(20);
        
    else
        disp(['Result already exists, calculation skipped... ' num2str(foldset) ' of ' num2str(Nfold)])
    end
    
    save([logpth '/D' num2str(tag) num2str(subtag) '_' num2str(foldset) 'of' num2str(Nfold) '.mat'],'foldset');
    waitbar(foldset/Nfold,wh);
    
end
close(wh);

%%
d=sprintf(['\nFinished ' Exp_str '-train on Data' num2str(tag) subtag ': loaded rcdt space data; performed PLDA on training data; saved PLDA directions, training data projections on the directions, and visualization axis.\n']); disp(d);
