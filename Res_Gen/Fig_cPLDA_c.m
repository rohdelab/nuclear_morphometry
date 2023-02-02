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
SD_spread=3;
which_axis=which_axisC;

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold which_axis reg_str reg_str_both
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


%                         tag=502;
%                         subtag='_test1';
%                         DISEASE_DETAIL={};
%                         DISEASE_SHORT={
%                             'FHB','FNH','HCA'
%                             };


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
            nme=['Dataset_' num2str(tr_tag)];
            %


            tr_ote={'1a','1b','1c','1d'}; tr_otnme={'liver','thyroid','mesothelioma','melanoma'};
            tr_Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};

        otherwise
            disp('Terminating...')
    end

    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/MAIN'];
    Exp_str = 'Main'; Exp_strF = '/MAIN/';
    cd(p0);

    %%
    wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
    for foldset=1:1%Nfold
        if reg_str_both==1
            vnm1=[Exp_str '_data' num2str(tr_tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=[Exp_str '_data' num2str(tr_tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=[Exp_str '_data' num2str(tr_tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        else
            vnm1=[Exp_str '_data' num2str(tr_tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=[Exp_str '_data' num2str(tr_tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=[Exp_str '_data' num2str(tr_tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        end

        vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];

        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);

        if ISthere1==0|ISthere2==0|ISthere3==0
            disp('ERROR!! Result doesnot exist...'); break;
        else
            disp([Exp_str ' testC on Data-' num2str(tr_tag) '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split

            load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);

            inp_te=[pp '/DATA/data' num2str(tag) subtag '/lotp'];

            load([inp_te '/Lotp_' reg_str]);
            data_te=u; label_te=label;

            size_multiplier=1.2;
            [fh_all1,fh_all2]=proj_n_plot_PLDAl_testC_ResGen(data_te,label_te,PLDA_directions,viz_plda,which_axis,labelnames,size_multiplier);


            Fgflnm=['Fig_' 'cPLDA_c' '/'];
            if exist(Fgflnm)
            else
                mkdir(Fgflnm)
            end

            for a=1:length(fh_all1)
                f01=fh_all1{a};
                if a==length(fh_all1)
                    %                     saveas(f01,[Fgflnm 'c1_axALL_tag_' num2str(tag) '.svg']);
                    saveas(f01,[Fgflnm 'c1_axALL_tag_' num2str(tag) '.png']);
                else
                    %                     saveas(f01,[Fgflnm 'c1_ax' num2str(which_axis(a)) 'tag_' num2str(tag) '.svg']);
                    saveas(f01,[Fgflnm 'c1_ax' num2str(which_axis(a)) 'tag_' num2str(tag) '.png']);
                end
                pause(1); close(f01);
            end
            for a=1:length(fh_all2)
                f01=fh_all2{a};

                %                 saveas(f01,[Fgflnm 'c2_ax' num2str(which_axis(a)) 'tag_' num2str(tag) '.svg']);
                saveas(f01,[Fgflnm 'c2_ax' num2str(which_axis(a)) 'tag_' num2str(tag) '.png']);

                pause(1); close(f01);
            end

            pause(1);

            disp(' ');
        end
        waitbar(foldset/Nfold,wh);
    end

    close(wh); close all;
    d=sprintf(['\nFinished ' Exp_str '-testC on Data' num2str(tag) subtag '__: projected test data on PLDA directions in LOTP space; saved PLDA projections, and 1D 2D plots.\n\n']); disp(d);
end
%%
