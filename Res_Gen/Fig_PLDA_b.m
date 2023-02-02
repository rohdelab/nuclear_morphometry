clc
clear all
close all
warning('off','all')

%% PARAMETER LOAD - selected in main.m
TagVec = 1; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=10; Nfold=2; Ndir01=5;

reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'
reg_str_both=1; % use both 'TOF' and 'TSOF' if set = 1

which_axisA=[1 2]; which_axisC=[1:4];  which_axisD=1;
pca_or_plda_in_01=0; % 1 = PCA, 0 = PLDA

which_axis=which_axisA; Ndir=Ndir01;

addpath('../codeM/')

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold which_axis reg_str reg_str_both pca_or_plda_in_01
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
        otherwise
            disp('Terminating...')
    end

    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/MAIN'];
    Exp_str = 'Main'; Exp_strF = '/MAIN/';
    cd(p0);

    if pca_or_plda_in_01==1
        str_save_dir='_PCA';
    else
        str_save_dir='_PLDA';
    end

    %%
    for foldset=1:1 % Nfold
        if reg_str_both==1
            vnm1=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        else
            vnm1=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=['main01' str_save_dir(2:end) Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        end

        vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];

        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);

        if ISthere1==0|ISthere2==0|ISthere3==0
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

                load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);

                inp_te=[pp '/DATA/data' ote{a} '/lotp'];
                inp_patient_label_te=[pp '/DATA/data' num2str(tag) subtag{a}];
                load([inp_te '/Lotp_' reg_str]); load([inp_patient_label_te '/patient_label' num2str(tag) subtag{a}]);
                data_te=u; label_te=label; patient_label_te{a}=label_patient;

                [xx_tr,lab_tr,xx_te,lab_te]=fn_trte_split2(data_te,label_te,ind,fold_tr,fold_te);
                [~,pat_lab_te]=fn_trte_split_patient_label(patient_label_te{a},ind,fold_tr,fold_te);

                size_multiplier=1;
                [projPLDA_te,f01,f02,~]=proj_n_plot_PLDAl_ResGen(xx_tr,lab_tr,xx_te,lab_te,PLDA_directions{a},viz_plda{a},which_axis,size_multiplier,Dnms{a});

                disp([Dnms{a} ' data:-']);
                [H_tmp,P_tmp] = mult_var_ttest2(lab_te,projPLDA_te);
                H_big(a)=H_tmp; P_big(a)=P_tmp;
                if H_tmp==0
                    disp(['Multivariate t-test; NOT statistically significant (!!); p-value (%): ' num2str(P_tmp*100)])
                else
                    disp(['Multivariate t-test; Statistically significant (\/); p-value (%): ' num2str(P_tmp*100)])
                end

                Fgflnm=['Fig_' 'PLDA' '/'];
                if exist(Fgflnm)
                else
                    mkdir(Fgflnm)
                end
                saveas(f01,[Fgflnm Dnms{a} '1.svg']);
                saveas(f02,[Fgflnm Dnms{a} '2.svg']);

                saveas(f01,[Fgflnm Dnms{a} '1.png']);
                saveas(f02,[Fgflnm Dnms{a} '2.png']);

                pause(1); close all;
            end
            disp(' ');
        end
    end

    close all;
end
%%
