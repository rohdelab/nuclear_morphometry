clc
clear all
close all
warning('off','all')

%% PARAMETER LOAD - selected in main.m
load(['../DATA/METADATA/nfe_params']); % Ndir, Nfold, reg_str, which_axisA, which_axisC, which_axisD, TagVec
which_axis=which_axisA; Ndir=Ndir01;

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold which_axis reg_str pca_or_plda_in_01
    %% SELECT DATA
    switch i_tag
        case 0
            tag=0;
            subtag{1}='y'; % toy data - 3
            subtag{2}='z'; % toy data - 2
            subtag{3}=[]; % toy data % subtag is a character like 'a', 'b', etc.

            nme=['Dataset_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
            Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};
        case 1
            tag=1;
            subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

            nme=['Dataset_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
            Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};
        case 2
            tag=2;
            subtag{1}='d';
            subtag{2}='d';

            nme=['Dataset_' num2str(tag)];
            ote={'2d','2d'}; otnme={'prostate','prostate'};
            Dnms={'Prostate','Prostate'};

        case 3
            tag=3;
            subtag{1}='a';
            subtag{2}='a';

            nme=['Dataset_' num2str(tag)];
            ote={'3a','3a'}; otnme={'meso_can_reac','meso_can_reac'};
            Dnms={'Meso_can_reac','Meso_can_reac'};

        otherwise
            disp('Terminating...')
    end

    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/MAIN_CMP_NFE'];
    Exp_str = 'Main_cmp_nfe'; Exp_strF = '/MAIN_CMP_NFE/';
    cd(p0);

    if pca_or_plda_in_01==1
        str_save_dir='_PCA';
    else
        str_save_dir='_PLDA';
    end

    %%
    wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
    for foldset=1:Nfold
        vnm1=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
        vnm2=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
        vnm3=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];


        vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];

        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);

        if ISthere1==0|ISthere2==0
            disp('ERROR!! Result doesnot exist...'); break;
        else
            disp([Exp_str ' test01A on Data-' num2str(tag) '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split

            H_big=[]; P_big=[];
            for a=1:length(ote)
                vnm4=['main01' Exp_str '_data' num2str(tag) 'ab' ote{a} '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
                vnm4pnm=[mdpth Exp_strF vnm4];

                indnm=['run5_indsplit_data' ote{a} '_fold' num2str(Nfold)];
                indpnm=[mdpth '/' indnm];
                load(indpnm)

                load(vnm1pnm); load(vnm2pnm); % load(vnm3pnm);

                inp_te=[pp '/DATA/data' ote{a} '/nfe'];
                inp_patient_label_te=[pp '/DATA/data' num2str(tag) subtag{a}];
                load([inp_te '/Nfe_' reg_str]); load([inp_patient_label_te '/patient_label' num2str(tag) subtag{a}]);
                data_te=u; label_te=label; patient_label_te{a}=label_patient;

                [xx_tr,lab_tr,xx_te,lab_te]=fn_trte_split2(data_te,label_te,ind,fold_tr,fold_te);
                [~,pat_lab_te]=fn_trte_split_patient_label(patient_label_te{a},ind,fold_tr,fold_te);


                [projPLDA_te]=proj_n_plot_PLDAl_cmp_nfe(xx_tr,lab_tr,xx_te,lab_te,PLDA_directions{a},which_axis);

                ESC=char(27); [ESC '[35m']; disp([Dnms{a} ' data:-']);
                [H_tmp,P_tmp] = mult_var_ttest2(lab_te,projPLDA_te);
                H_big(a)=H_tmp; P_big(a)=P_tmp;
                if H_tmp==0
                    disp(['Multivariate t-test; NOT statistically significant (!!); p-value (%): ' num2str(P_tmp*100)])
                else
                    disp(['Multivariate t-test; Statistically significant (\/); p-value (%): ' num2str(P_tmp*100)])
                end

                [MM_all,MM_label_all,HH_all] = patient_hist_mean_feature(projPLDA_te,lab_te,pat_lab_te,1);

                tmp_pth=[mdpth Exp_strF];
                tmp_nm1=[str_save_dir(2:end) 'Feat1_test' '_data' ote{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
                tmp_nm2=[str_save_dir(2:end) 'Feat2_test' '_data' ote{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
                tmp_nm3=[str_save_dir(2:end) 'Feat3_test' '_data' ote{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];


                feat=MM_all';
                y_pat=MM_label_all; y_pat(y_pat>2)=2;
                save([tmp_pth tmp_nm1],'feat','y_pat');

                feat=projPLDA_te';
                y_cell=lab_te; y_cell(y_cell>2)=2;
                y_pat=pat_lab_te;
                save([tmp_pth tmp_nm2],'feat','y_cell','y_pat');

                feat=HH_all';
                y_pat=MM_label_all; y_pat(y_pat>2)=2;
                save([tmp_pth tmp_nm3],'feat','y_pat');


                save(vnm4pnm,'projPLDA_te','-v7.3');


                if exist([respth '/' nme '/test01A' str_save_dir '/'])
                else
                    mkdir([respth '/' nme '/test01A' str_save_dir '/'])
                end

            end
            disp(' ');
        end
        waitbar(foldset/Nfold,wh);
    end

    close(wh); close all;
    d=sprintf(['\nFinished ' Exp_str '-test01A on Data' num2str(tag) '__: projected test data on PLDA directions in LOTP space; saved PLDA projections, and 1D 2D plots.\n\n']); disp(d);
end
%%
