clc
clear all
close all
warning('off','all')

%% PARAMuETER SELECT
TagVec = 1; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=5; Nfold=2;
which_axis=[1 2]; % entries must be <= Ndir, for testA only
reg_str='TOF'; % 'TOF' 'TSOF'

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold which_axis reg_str
    %% SELECT DATA
    switch i_tag
        case 1
            %             tag=501;
            %             subtag='_test1';
            %             DISEASE_DETAIL={
            %                 'Breast invasive carcinoma','Breast invasive carcinoma','Breast invasive carcinoma','Breast invasive carcinoma','Breast invasive carcinoma',...
            %                 'Breast invasive carcinoma','Kidney renal clear cell carcinoma','Kidney renal papillary cell carcinoma','Kidney renal papillary cell carcinoma','Kidney renal papillary cell carcinoma',...
            %                 'Kidney renal clear cell carcinoma','Kidney renal clear cell carcinoma','Liver/Lung squamous cell carcinoma','Liver/Lung adenocarcinoma','Liver/Lung adenocarcinoma',...
            %                 'Liver/Lung adenocarcinoma','Liver/Lung squamous cell carcinoma','Liver/Lung squamous cell carcinoma','Prostate adenocarcinoma','Prostate adenocarcinoma',...
            %                 'Prostate adenocarcinoma','Prostate adenocarcinoma','Prostate adenocarcinoma','Prostate adenocarcinoma','Bladder Urothelial Carcinoma',...
            %                 'Bladder Urothelial Carcinoma','Colon adenocarcinoma','Colon adenocarcinoma','Stomach adenocarcinoma','Stomach adenocarcinoma'
            %                 };
            %             DISEASE_SHORT={
            %                 'Breast','Breast','Breast','Breast','Breast',...
            %                 'Breast','Kidney','Kidney','Kidney','Kidney',...
            %                 'Kidney','Kidney','Liver/Lung','Liver/Lung','Liver/Lung',...
            %                 'Liver/Lung','Liver/Lung','Liver/Lung','Prostate','Prostate',...
            %                 'Prostate','Prostate','Prostate','Prostate','Bladder',...
            %                 'Bladder','Colon','Colon','Stomach','Stomach',...
            %                 };
            
            
            %             tag=501;
            %             subtag='_test2';
            %             DISEASE_DETAIL={};
            %             DISEASE_SHORT={
            %                 'Breast','Kidney','Liver','Prostate','Bladder',...
            %                 'Colon','Stomach'
            %                 };
            
            
            %             tag=502;
            %             subtag='_test1';
            %             DISEASE_DETAIL={};
            %             DISEASE_SHORT={
            %                 'FHB','FNH','HCA'
            %                 };
            
            
            tag=503;
            subtag='_test1';
            DISEASE_DETAIL={};
            DISEASE_SHORT={
                'FA','FC','FVPC','NG'
                };
            
            
            
            % % % % % % % % % % % % %
            
            labelnames=DISEASE_SHORT;
            
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            tr_tag=1;
            tr_subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            tr_subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            tr_subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            tr_subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
            
            %
            %             nme=['D_composite_te_' num2str(tag) subtag '_tr_' num2str(tr_tag)];
            nme=['D_composite_' num2str(tr_tag)];
            %
            tr_ote={'1a','1b','1c','1d'}; tr_otnme={'liver','thyroid','mesothelioma','melanoma'};
            tr_Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};
            
        otherwise
            disp('Terminating...')
    end
    
    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/EXP05_ex2'];
    Exp_str = 'Exp05_ex2'; Exp_strF = '/EXP05_ex2/';
    cd(p0);
    
    %%
    wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
    for foldset=1:Nfold
        vnm1=[Exp_str '_data' num2str(tr_tag) tr_subtag{1} tr_subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
        vnm2=[Exp_str '_data' num2str(tr_tag) tr_subtag{1} tr_subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
        vnm3=[Exp_str '_data' num2str(tr_tag) tr_subtag{1} tr_subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        
        vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];
        
        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
        
        if ISthere1==0|ISthere2==0|ISthere3==0
            disp('ERROR!! Result doesnot exist...'); break;
        else
            disp([Exp_str ' testC on Data-' num2str(tr_tag) tr_subtag{1} tr_subtag{2} tr_subtag{3} tr_subtag{4} '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split
            
            load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
            
            inp_te=[pp '/DATA/data' num2str(tag) subtag '/lotp'];
            
            load([inp_te '/Lotp_' reg_str]);
            data_te=u; label_te=label;
            
            [fh_all]=proj_n_plot_PLDAl_testC(data_te,label_te,PLDA_directions,viz_plda,which_axis,labelnames);
            
            %
            
            %             if exist([respth '/' nme '/testC/'])
            %             else
            %                 mkdir([respth '/' nme '/testC/'])
            %             end
            %
            %             for a=1:length(fh_all)
            %                 f01=fh_all{a};
            %                 saveas(f01,[respth '/' nme '/testC/' num2str(tag) subtag '_ax' num2str(which_axis(a)) '_set' num2str(foldset) '_' reg_str '.png']);
            %                 pause(10);
            %             end
            
            if exist([respth '/' nme '/testC/' num2str(tag) subtag '/'])
            else
                mkdir([respth '/' nme '/testC/' num2str(tag) subtag '/'])
            end
            
            for a=1:length(fh_all)
                f01=fh_all{a};
                saveas(f01,[respth '/' nme '/testC/' num2str(tag) subtag '/ax' num2str(which_axis(a)) '_set' num2str(foldset) '_' reg_str '.png']);
                pause(10);
            end
            
            
            
            %
            close all;
            
            disp(' ');
        end
        waitbar(foldset/Nfold,wh);
    end
    
    close(wh); close all;
    d=sprintf(['\nFinished ' Exp_str '-testC on Data' num2str(tag) subtag '__: projected test data on PLDA directions in LOTP space; saved PLDA projections, and 1D 2D plots.\n\n']); disp(d);
end
%%
