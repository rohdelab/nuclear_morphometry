clc
clear all
close all
warning('off','all')

%% PARAMETER LOAD - selected in main.m
load(['../DATA/METADATA/nfe_params']); % Ndir, Nfold, reg_str, which_axisA, which_axisC, which_axisD, TagVec, reg_str_both
Ndir=Ndir01;

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold reg_str reg_str_both pca_or_plda_in_01
    %% SELECT DATA
    switch i_tag
        case 0
            tag=0;
            subtag{1}='y'; % toy data - 3
            subtag{2}='z'; % toy data - 2
            subtag{3}=[]; % toy data % subtag is a character like 'a', 'b', etc.
        case 1
            tag=1;
            subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
        case 2
            tag=2;
            subtag{1}='d';
            subtag{2}='d';
        case 3
            tag=3;
            subtag{1}='a';
            subtag{2}='a';
        otherwise
            disp('Terminating...')
    end

    %% LOAD DATA & PRELIMINARIES
    p0=pwd; cd ..
    mdpth=[pwd '/DATA/METADATA'];

    for a=1:length(subtag)
        inp_tr=[pwd '/DATA/data' num2str(tag) subtag{a} '/nfe'];
        inp_patient_label_tr=[pwd '/DATA/data' num2str(tag) subtag{a}];

        load([inp_tr '/Nfe_' 'TOF']); load([inp_patient_label_tr '/patient_label' num2str(tag) subtag{a}]);
        data_tr1{a}=u; label_tr1{a}=label; patient_label_tr1{a}=label_patient;

        load([inp_tr '/Nfe_' 'TOF']); load([inp_patient_label_tr '/patient_label' num2str(tag) subtag{a}]); % not a mistake - TOF
        data_tr2{a}=u; label_tr2{a}=label; patient_label_tr2{a}=label_patient;
    end

    Exp_str = 'Main_cmp_nfe'; Exp_strF = '/MAIN_CMP_NFE/';
    cd(p0);

    if pca_or_plda_in_01==1
        str_save_dir='_PCA';
    else
        str_save_dir='_PLDA';
    end

    %% CALCULATION
    wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
    for foldset=1:Nfold
        tic

        vnm1=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
        vnm2=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
        vnm3=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];

        if exist([mdpth Exp_strF])
        else
            mkdir([mdpth Exp_strF])
        end

        vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];

        if pca_or_plda_in_01==1
            disp(['PCA-train on Data-' num2str(tag) '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
        else
            disp(['PLDA-train on Data-' num2str(tag) '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
        end

        fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split

        xxT1=[]; labT1=[];
        for a=1:length(subtag)
            indnm=['run5_indsplit_data' num2str(tag) subtag{a} '_fold' num2str(Nfold)];
            indpnm=[mdpth '/' indnm];
            load(indpnm)
            [xx_tr_s,lab_tr_s,~,~]=fn_trte_split2(data_tr1{a},label_tr1{a},ind,fold_tr,fold_te);
            [pat_lab_tr_s,~]=fn_trte_split_patient_label(patient_label_tr1{a},ind,fold_tr,fold_te);
            xx_tr1{a}=xx_tr_s; lab_tr1{a}=lab_tr_s; pat_lab_tr1{a}=pat_lab_tr_s;
            xxT1=[xxT1 xx_tr_s]; ltmp=lab_tr_s; ltmp(ltmp>2)=2; labT1=[labT1;ltmp(:)];
        end

        xxT2=[]; labT2=[];
        for a=1:length(subtag)
            indnm=['run5_indsplit_data' num2str(tag) subtag{a} '_fold' num2str(Nfold)];
            indpnm=[mdpth '/' indnm];
            load(indpnm)
            [xx_tr_s,lab_tr_s,~,~]=fn_trte_split2(data_tr2{a},label_tr2{a},ind,fold_tr,fold_te);
            [pat_lab_tr_s,~]=fn_trte_split_patient_label(patient_label_tr2{a},ind,fold_tr,fold_te);
            xx_tr2{a}=xx_tr_s; lab_tr2{a}=lab_tr_s; pat_lab_tr2{a}=pat_lab_tr_s;
            xxT2=[xxT2 xx_tr_s]; ltmp=lab_tr_s; ltmp(ltmp>2)=2; labT2=[labT2;ltmp(:)];
        end

        if pca_or_plda_in_01==1
            [dir1,proj1]=scPCA_ModesL_ex1v2_both_cmp_nfe(xx_tr1,lab_tr1,xx_tr1,lab_tr1,Ndir,xxT1,labT1,xxT2,labT2);
        else
            [dir1,proj1]=scPLDA_ModesL_ex1v2_both_cmp_nfe(xx_tr1,lab_tr1,xx_tr1,lab_tr1,Ndir,xxT1,labT1,xxT2,labT2);
        end

        for a=1:length(proj1)
            PLDA_directions{a}=[dir1{a}];
            projPLDA_tr{a}=[proj1{a}];
        end

        lab_tr=lab_tr1; pat_lab_tr=pat_lab_tr1;

        [MM_all,MM_label_all,HH_all] = patient_hist_mean_feature(projPLDA_tr,lab_tr,pat_lab_tr,length(projPLDA_tr));

        for a=1:length(subtag)
            tmp_pth=[mdpth Exp_strF];
            tmp_nm1=[str_save_dir(2:end) 'Feat1_train' '_data' num2str(tag) subtag{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
            tmp_nm2=[str_save_dir(2:end) 'Feat2_train' '_data' num2str(tag) subtag{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
            tmp_nm3=[str_save_dir(2:end) 'Feat3_train' '_data' num2str(tag) subtag{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];

            feat=MM_all{a}';
            y_pat=MM_label_all{a}; y_pat(y_pat>2)=2;
            save([tmp_pth tmp_nm1],'feat','y_pat');

            feat=projPLDA_tr{a}';
            y_cell=lab_tr{a}; y_cell(y_cell>2)=2;
            y_pat=pat_lab_tr{a};
            save([tmp_pth tmp_nm2],'feat','y_cell','y_pat');

            feat=HH_all{a}';
            y_pat=MM_label_all{a}; y_pat(y_pat>2)=2;
            save([tmp_pth tmp_nm3],'feat','y_pat');
        end

        save(vnm1pnm,'PLDA_directions','-v7.3');
        save(vnm2pnm,'projPLDA_tr','lab_tr','-v7.3');
        pause(1);

        waitbar(foldset/Nfold,wh);
        toc
        disp(' ');
    end

    close(wh);
    d=sprintf(['\nFinished ' Exp_str '-train 01 on Data' num2str(tag) 'ab' '__: loaded lotp space data; performed PLDA on training data; saved PLDA directions, training data projections on the directions, and visualization axis.\n\n']); disp(d);
end

