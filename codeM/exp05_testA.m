clc
clear all
close all
warning('off','all')

%% PARAMuETER SELECT
TagVec = 1; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=50; Nfold=2;
which_axis=[1 2]; % entries must be <= Ndir, for testA only
reg_str='TOF'; % 'TOF' 'TSOF'

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold which_axis reg_str
    %% SELECT DATA
    switch i_tag
        case 0
            tag=0;
            subtag{1}='y'; % toy data - 3
            subtag{2}='z'; % toy data - 2
            subtag{3}=[]; % toy data % subtag is a character like 'a', 'b', etc.
            
            nme=['D_composite_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
            Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};
        case 1
            tag=1;
            subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
            
            nme=['D_composite_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
            Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};
        otherwise
            disp('Terminating...')
    end
    
    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/EXP05'];
    Exp_str = 'Exp05'; Exp_strF = '/EXP05/';
    cd(p0);
    
    %%
    wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
    for foldset=1:Nfold
        vnm1=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
        vnm2=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
        vnm3=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        vnm5=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_pat_hist_mean_class_tr'];
        
        vnm1pnm=[mdpth Exp_strF vnm1];
        vnm2pnm=[mdpth Exp_strF vnm2];
        vnm3pnm=[mdpth Exp_strF vnm3];
        vnm5pnm=[mdpth Exp_strF vnm5];
        
        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
        
        if ISthere1==0|ISthere2==0|ISthere3==0
            disp('ERROR!! Result doesnot exist...'); break;
        else
            disp([Exp_str ' testA on Data-' num2str(tag) subtag{1} subtag{2} subtag{3} subtag{4} '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split
            
            %             disp('Patient classification (hist mean feature based); Test accuracy:')
            H_big=[]; P_big=[];
            %
            inT=[pp '/DATA/data' num2str(tag) 'T' '/lotp'];
            load([inT '/Lotp_' reg_str]);
            dT=u; pT=ptcl_wght; lT=label;
            st=1;en=0;
            %
            for a=1:length(ote)
                vnm4=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} ote{a} '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
                vnm4pnm=[mdpth Exp_strF vnm4];
                
                indnm=['run5_indsplit_data' ote{a} '_fold' num2str(Nfold)];
                indpnm=[mdpth '/' indnm];
                load(indpnm)
                
                load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);

                % %
                inp_patient_label_te=[pp '/DATA/data' num2str(tag) subtag{a}];
                load([inp_patient_label_te '/patient_label' num2str(tag) subtag{a}]);
                patient_label_te{a}=label_patient;
                
                en=en+length(label_patient);
                data_te=dT(:,st:en); label_te=lT(st:en);
                st=1+en;
                % %
                
                
                [xx_tr,lab_tr,xx_te,lab_te]=fn_trte_split2(data_te,label_te,ind,fold_tr,fold_te);
                [~,pat_lab_te]=fn_trte_split_patient_label(patient_label_te{a},ind,fold_tr,fold_te);
                [projPLDA_te,f01,f02,~]=proj_n_plot_PLDAl(xx_tr,lab_tr,xx_te,lab_te,PLDA_directions,viz_plda,which_axis);
                
                ESC=char(27); [ESC '[35m']; disp([Dnms{a} ' data:-']);
                [H_tmp,P_tmp] = mult_var_ttest2(lab_te,projPLDA_te);
                H_big(a)=H_tmp; P_big(a)=P_tmp;
                if H_tmp==0
                    disp(['Multivariate t-test; '...
                        ESC '[33m' 'NOT statistically significant (!!); p-value (%): ' num2str(P_tmp*100) ESC '[35m'])
                else
                    disp([ESC '[32m' 'Multivariate t-test; '...
                        ESC '[33m' 'Statistically significant (\/); p-value (%): ' num2str(P_tmp*100) ESC '[35m'])
                end
                
                % % %
                
                %                 [MM_all,MM_label_all] = patient_hist_mean_feature(projPLDA_te,lab_te,pat_lab_te,1);
                [MM_all,MM_label_all,HH_all] = patient_hist_mean_feature(projPLDA_te,lab_te,pat_lab_te,1);
                
                tmp_pth=[mdpth Exp_strF];
                tmp_nm1=['Feat1_test' '_data' ote{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
                tmp_nm2=['Feat2_test' '_data' ote{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
                tmp_nm3=['Feat3_test' '_data' ote{a} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir)];
                
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
                
                
                % % %
                
                X=MM_all'; y_true=MM_label_all;
                load(vnm5pnm);
                y_pred = svmclassify(svm_trained{a},X,'showplot',true);
                phmc_acc_tr(foldset,a)=100*length(find(y_true(:)==y_pred(:)))/length(y_true);
                disp([ESC '[32m' 'Patient classification (hist mean feature based); '...
                    ESC '[33m' 'Test accuracy: ' num2str(phmc_acc_tr(foldset,a)) ' %' ESC '[35m']);
                
                save(vnm4pnm,'projPLDA_te','-v7.3');
                saveas(f01,[respth '/' nme '/testA/f1_' otnme{a} '_ax' num2str(which_axis(1)) '_set' num2str(foldset) '_' reg_str '.png']);
                saveas(f02,[respth '/' nme '/testA/f2_' otnme{a} '_ax' num2str(which_axis(1)) '_' num2str(which_axis(2)) '_set' num2str(foldset) '_' reg_str '.png']);
                pause(10); close all;
            end
            disp(' ');
        end
        waitbar(foldset/Nfold,wh);
    end
    
    disp(' ');
    phmc_acc_tr_mean=mean(phmc_acc_tr);
    phmc_acc_tr_std=std(phmc_acc_tr);
    disp(['mean accuracy (%): ' num2str(phmc_acc_tr_mean)]);
    disp(['std accuracy (%): ' num2str(phmc_acc_tr_std)]);
    
    close(wh); close all;
    d=sprintf(['\nFinished ' Exp_str '-testA on Data' num2str(tag) subtag{1} subtag{2} subtag{3} subtag{4} '__: projected test data on PLDA directions in LOTP space; saved PLDA projections, and 1D 2D plots.\n\n']); disp(d);
end
%%
