clc
clear all
close all
warning('off','all')

%% PARAMETER SELECT
TagVec = 1; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=5; Nfold=2;
force_run=1;% 1 for must run no matter what, 0 for may skip if already there.
reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold force_run reg_str
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
        otherwise
            disp('Terminating...')
    end
    
    %% LOAD DATA & PRELIMINARIES
    p0=pwd; cd ..
    mdpth=[pwd '/DATA/METADATA'];
    logpth=[pwd '/RESULTS/EXP06_ex1/log'];
    dT=[]; pT=[]; lT=[];
    for a=1:length(subtag)
        inp_tr=[pwd '/DATA/data' num2str(tag) subtag{a} '/lotp'];
        inp_patient_label_tr=[pwd '/DATA/data' num2str(tag) subtag{a}];
        load([inp_tr '/Lotp_' reg_str]); load([inp_patient_label_tr '/patient_label' num2str(tag) subtag{a}]);
        data_tr{a}=u; p_wt_tr{a}=ptcl_wght; label_tr{a}=label; patient_label_tr{a}=label_patient;
    end
    Exp_str = 'Exp06_ex1'; Exp_strF = '/EXP06_ex1/';
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
            
            ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
            
            if force_run|ISthere1==0|ISthere2==0|ISthere3==0
                disp([Exp_str '-train on Data-' num2str(tag) '_exclude - ' subtag{ii} '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
                
                fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split
                
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
                
                
                [PLDA_directions,projPLDA_tr,viz_plda]=cPLDA_ModesL_ex1(xx_tr,lab_tr,p_wt_tr,Ndir,xxT,pT,labT);
                
                
                save(vnm1pnm,'PLDA_directions','-v7.3');
                save(vnm2pnm,'projPLDA_tr','lab_tr','-v7.3');
                save(vnm3pnm,'viz_plda','-v7.3');
                pause(20);
            else
                disp(['Result already exists, calculation skipped... ' num2str(foldset) ' of ' num2str(Nfold)])
            end
            
            if exist([logpth '/'])
            else
                mkdir([logpth '/'])
            end
            save([logpth '/D' num2str(tag) '_' num2str(foldset) 'of' num2str(Nfold) '.mat'],'foldset');
            waitbar(foldset/Nfold,wh);
            toc
            disp(' ');
        end
        
        

        
        close(wh);
        
    end
    d=sprintf(['\nFinished ' Exp_str '-train on Data' num2str(tag) subtag{1} subtag{2} '_exclude - ' subtag{ii} '__: loaded lotp space data; performed PLDA on training data; saved PLDA directions, training data projections on the directions, and visualization axis.\n\n']); disp(d);
end

