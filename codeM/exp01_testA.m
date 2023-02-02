clc
clear all
close all

warning('off','all')

%%
tag=0; subtag=[]; nme='D_toy'; ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'}; % toy data % subtag is a character like 'a', 'b', etc.

% tag=1; subtag='a'; nme='D_liver'; ote={'1b','1c','1d'}; otnme={'thyroid','mesothelioma','melanoma'}; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; nme='D_thyroid'; ote={'1a','1c','1d'}; otnme={'liver','mesothelioma','melanoma'}; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; nme='D_mesothelioma'; ote={'1a','1b','1d'}; otnme={'liver','thyroid','melanoma'}; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; nme='D_melanoma'; ote={'1a','1b','1c'}; otnme={'liver','thyroid','mesothelioma'}; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

%%
Ndir=20; Nfold=10;
which_axis=1; % must be <= Ndir - 1
reg_str='TOF';

%%
p0=pwd; cd ..; pp=pwd;
mdpth=[pwd '/DATA/METADATA'];
respth=[pwd '/RESULTS/EXP01'];
Exp_str = 'Exp01'; Exp_strF = '/EXP01/'; 
cd(p0);

%%
wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
for foldset=1:Nfold
    vnm1=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
    vnm2=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
    vnm3=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
    
    
    vnm1pnm=[mdpth Exp_strF vnm1];
    vnm2pnm=[mdpth Exp_strF vnm2];
    vnm3pnm=[mdpth Exp_strF vnm3];

    
    ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
    
    
    if ISthere1==0|ISthere2==0|ISthere3==0
        disp('ERROR!! Result doesnot exist...'); break;
    else
        disp([Exp_str ' testA on Data ' num2str(tag) subtag ': Calculating and saving... ' num2str(foldset) ' of ' num2str(Nfold)]);
        fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 10-fold split
        
        vnm4=[Exp_str '_data' num2str(tag) subtag num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
        vnm4pnm=[mdpth Exp_strF vnm4];
        
        indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
        indpnm=[mdpth '/' indnm];
        
        load(indpnm)
        load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
        
        inp_te=[pp '/DATA/data' num2str(tag) subtag '/rcdt'];
        load([inp_te '/Rcdt']);
        data_te=u; label_te=label;
        
        [~,~,xx_te,lab_te]=fn_trte_split(data_te,label_te,ind,fold_tr,fold_te);
        [projPLDA_te,f01,f02,~]=proj_n_plot_PLDAr(xx_te,lab_te,PLDA_directions,viz_plda,which_axis);
        
        save(vnm4pnm,'projPLDA_te','-v7.3');
        saveas(f01,[respth '/' nme '/testA/f1_' nme '_ax' num2str(which_axis) '_set' num2str(foldset) '.png']);
        saveas(f02,[respth '/' nme '/testA/f2_' nme '_ax' num2str(which_axis) '_set' num2str(foldset) '.png']);
        pause(10); close all;
        
        for a=1:length(ote)
            vnm4=[Exp_str '_data' num2str(tag) subtag ote{a} '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
            vnm4pnm=[mdpth Exp_strF vnm4];
            
            indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
            indpnm=[mdpth '/' indnm];
            
            load(indpnm)
            load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
            
            inp_te=[pp '/DATA/data' ote{a} '/rcdt'];
            load([inp_te '/Rcdt_' reg_str]);
            data_te=u; label_te=label;
            xx_te=data_te;  lab_te=label_te;
            [projPLDA_te,f01,f02,~]=proj_n_plot_PLDAr(xx_te,lab_te,PLDA_directions,viz_plda,which_axis);
            
            save(vnm4pnm,'projPLDA_te','-v7.3');
            saveas(f01,[respth '/' nme '/testA/f1_' otnme{a} '_ax' num2str(which_axis) '_set' num2str(foldset) '_' reg_str '.png']);
            saveas(f02,[respth '/' nme '/testA/f2_' otnme{a} '_ax' num2str(which_axis) '_set' num2str(foldset) '_' reg_str '.png']);
            pause(10); close all;
        end
        
    end
    
    waitbar(foldset/Nfold,wh);
end
close(wh); close all;

%%
d=sprintf(['\nFinished ' Exp_str '-testA on Data' num2str(tag) subtag ': projected test data on PLDA directions in RCDT space; saved PLDA projections, and 1D 2D plots.\n']); disp(d);
