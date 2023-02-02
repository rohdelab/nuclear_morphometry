clc
clear all
close all
warning('off','all')

%% PARAMETER LOAD - selected in main.m
load(['../DATA/METADATA/params']); % Ndir, Nfold, reg_str, which_axisA, which_axisC, which_axisD, TagVec, reg_str_both

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold reg_str reg_str_both
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
        
        otherwise
            disp('Terminating...')
    end

    %% LOAD DATA & PRELIMINARIES
    p0=pwd; cd ..
    mdpth=[pwd '/DATA/METADATA'];

    if reg_str_both==1
        reg_str='BOTH';
        for a=1:length(subtag)
            inp_tr=[pwd '/DATA/data' num2str(tag) subtag{a} '/lotp'];
            inp_patient_label_tr=[pwd '/DATA/data' num2str(tag) subtag{a}];

            load([inp_tr '/Lotp_' 'TOF']); load([inp_patient_label_tr '/patient_label' num2str(tag) subtag{a}]);
            data_tr1{a}=u; p_wt_tr1{a}=ptcl_wght; label_tr1{a}=label; patient_label_tr1{a}=label_patient;

            load([inp_tr '/Lotp_' 'TSOF']); load([inp_patient_label_tr '/patient_label' num2str(tag) subtag{a}]);
            data_tr2{a}=u; p_wt_tr2{a}=ptcl_wght; label_tr2{a}=label; patient_label_tr2{a}=label_patient;
        end
    else
        for a=1:length(subtag)
            inp_tr=[pwd '/DATA/data' num2str(tag) subtag{a} '/lotp'];
            inp_patient_label_tr=[pwd '/DATA/data' num2str(tag) subtag{a}];
            load([inp_tr '/Lotp_' reg_str]); load([inp_patient_label_tr '/patient_label' num2str(tag) subtag{a}]);
            data_tr{a}=u; p_wt_tr{a}=ptcl_wght; label_tr{a}=label; patient_label_tr{a}=label_patient;
        end
    end

    Exp_str = 'Main'; Exp_strF = '/MAIN/';
    cd(p0);

    %% CALCULATION
    for ii=1:length(subtag)
        wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
        for foldset=1:Nfold
            tic

            vnm1=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_excld' subtag{ii} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_excld' subtag{ii} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_excld' subtag{ii} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];

            if exist([mdpth Exp_strF])
            else
                mkdir([mdpth Exp_strF])
            end

            vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];

            disp(['Joint-PLDA-train on Data-' num2str(tag) '_exclude - ' subtag{ii} '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');

            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split

            if reg_str_both==1
                xxT1=[]; labT1=[]; pT1=[];
                cnt=0;
                for a=1:length(subtag)
                    indnm=['run5_indsplit_data' num2str(tag) subtag{a} '_fold' num2str(Nfold)];
                    indpnm=[mdpth '/' indnm];
                    load(indpnm)
                    if a==ii
                    else
                        cnt=cnt+1;
                        [xx_tr_s,lab_tr_s,~,~]=fn_trte_split2(data_tr1{a},label_tr1{a},ind,fold_tr,fold_te);
                        xx_tr1{cnt}=xx_tr_s; lab_tr1{cnt}=lab_tr_s;
                        xxT1=[xxT1 xx_tr_s]; ltmp=lab_tr_s; ltmp(ltmp>2)=2; labT1=[labT1;ltmp(:)]; pT1=p_wt_tr1{1};
                    end
                end

                xxT2=[]; labT2=[]; pT2=[];
                cnt=0;
                for a=1:length(subtag)
                    indnm=['run5_indsplit_data' num2str(tag) subtag{a} '_fold' num2str(Nfold)];
                    indpnm=[mdpth '/' indnm];
                    load(indpnm)
                    if a==ii
                    else
                        cnt=cnt+1;
                        [xx_tr_s,lab_tr_s,~,~]=fn_trte_split2(data_tr2{a},label_tr2{a},ind,fold_tr,fold_te);
                        xx_tr2{cnt}=xx_tr_s; lab_tr2{cnt}=lab_tr_s;
                        xxT2=[xxT2 xx_tr_s]; ltmp=lab_tr_s; ltmp(ltmp>2)=2; labT2=[labT2;ltmp(:)]; pT2=p_wt_tr2{1};
                    end
                end

                [dir1,proj1,viz1]=cPLDA_ModesL_ex1v2_both(xx_tr1,lab_tr1,p_wt_tr1,xx_tr1,lab_tr1,p_wt_tr1,Ndir,xxT1,pT1,labT1,xxT2,pT2,labT2);
                [dir2,proj2,viz2]=cPLDA_ModesL_ex1v2_both(xx_tr2,lab_tr2,p_wt_tr2,xx_tr1,lab_tr1,p_wt_tr1,Ndir,xxT2,pT2,labT2,xxT2,pT2,labT2);

                PLDA_directions=[dir1 dir2];
                for a=1:length(proj1)
                    projPLDA_tr{a}=[proj1{a};proj2{a}];
                end
                viz_plda=cat(3,viz1,viz2);

                lab_tr=lab_tr1;
            else
                xxT=[]; labT=[]; pT=[];
                cnt=0;
                for a=1:length(subtag)
                    indnm=['run5_indsplit_data' num2str(tag) subtag{a} '_fold' num2str(Nfold)];
                    indpnm=[mdpth '/' indnm];
                    load(indpnm)

                    if a==ii
                    else
                        cnt=cnt+1;
                        [xx_tr_s,lab_tr_s,~,~]=fn_trte_split2(data_tr{a},label_tr{a},ind,fold_tr,fold_te);
                        xx_tr{cnt}=xx_tr_s; lab_tr{cnt}=lab_tr_s;
                        xxT=[xxT xx_tr_s]; ltmp=lab_tr_s; ltmp(ltmp>2)=2; labT=[labT;ltmp(:)]; pT=p_wt_tr{1};
                    end
                end
                %                 [PLDA_directions,projPLDA_tr,viz_plda]=cPLDA_ModesL(xx_tr,lab_tr,p_wt_tr,Ndir,xxT,pT,labT);
                [PLDA_directions,projPLDA_tr,viz_plda]=cPLDA_ModesL_ex1v2(xx_tr,lab_tr,p_wt_tr,Ndir,xxT,pT,labT);
            end

            save(vnm1pnm,'PLDA_directions','-v7.3');
            save(vnm2pnm,'projPLDA_tr','lab_tr','-v7.3');
            save(vnm3pnm,'viz_plda','-v7.3');
            pause(1);

            waitbar(foldset/Nfold,wh);
            toc
            disp(' ');
        end

        close(wh);

    end
    d=sprintf(['\nFinished ' Exp_str '-train on Data' num2str(tag) subtag{1} subtag{2} '_exclude - ' subtag{ii} '__: loaded lotp space data; performed PLDA on training data; saved PLDA directions, training data projections on the directions, and visualization axis.\n\n']); disp(d);
end

