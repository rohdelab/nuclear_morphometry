clc
clear all
close all
warning('off','all')

%% PARAMuETER SELECT
TagVec = 1; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=5; Nfold=2;
which_axis=1; %[1:4]; % entries must be <= Ndir, for testA only
reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold which_axis reg_str
    %% SELECT DATA
    switch i_tag
        case 1
            tr_tag=1;
            tr_subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            tr_subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            tr_subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            tr_subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
            
            
            tr_ote={'1a','1b','1c','1d'}; tr_otnme={'liver','thyroid','mesothelioma','melanoma'};
            tr_Dnms={'Liver','Thyroid','Mesothelioma','Melanoma'};
            
        otherwise
            disp('Terminating...')
    end
    
    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/EXP06_ex1'];
    Exp_str = 'Exp06_ex1'; Exp_strF = '/EXP06_ex1/';
    cd(p0);
    
    %%
    for ii=1:length(tr_subtag)
        tag=tr_tag;
        subtag=tr_subtag{ii};
        DISEASE_DETAIL={};
        DISEASE_SHORT={
            'Less malignant','More malignant'
            };
        labelnames=DISEASE_SHORT;
        nme=['D_composite_' num2str(tr_tag)];
        
        wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
        for foldset=1:Nfold
            vnm1=[Exp_str '_data' num2str(tr_tag) tr_subtag{1} tr_subtag{2} '_excld' tr_subtag{ii} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=[Exp_str '_data' num2str(tr_tag) tr_subtag{1} tr_subtag{2} '_excld' tr_subtag{ii} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=[Exp_str '_data' num2str(tr_tag) tr_subtag{1} tr_subtag{2} '_excld' tr_subtag{ii} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
            
            vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];
            
            ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
            
            if ISthere1==0|ISthere2==0|ISthere3==0
                disp('ERROR!! Result doesnot exist...'); break;
            else
                disp([Exp_str ' testD on Data-' num2str(tr_tag) tr_subtag{1} tr_subtag{2} tr_subtag{3} tr_subtag{4} '_excld' tr_subtag{ii} '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
                fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split
                
                load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
                
                inp_te=[pp '/DATA/data' num2str(tag) subtag '/lotp'];
                
                load([inp_te '/Lotp_' reg_str]);
                data_te=u; label_te=label;
                
                disp([tr_Dnms{ii} ' data in testing; others in training ...']);
                [fh_all]=proj_n_plot_PLDAl_testD(data_te,label_te,PLDA_directions,viz_plda,which_axis,labelnames);

                debug_stop=4;
                
                if exist([respth '/' nme '/testD/' num2str(tag) subtag '/'])
                else
                    mkdir([respth '/' nme '/testD/' num2str(tag) subtag '/'])
                end
                
                for a=1:length(fh_all)
                    f01=fh_all{a};
                    saveas(f01,[respth '/' nme '/testD/' num2str(tag) subtag '/ax' num2str(which_axis(a)) '_set' num2str(foldset) '_' reg_str '.png']);
                    pause(10);
                end
                close all;
                
                disp(' ');
            end
            waitbar(foldset/Nfold,wh);
        end
        
        close(wh); pause(20); close all;
    end
    d=sprintf(['\nFinished ' Exp_str '-testD on Data' num2str(tag) subtag '_exclude - ' tr_subtag{ii} '__: projected test data on PLDA directions in LOTP space; saved PLDA projections, and 1D 2D plots.\n\n']); disp(d);
end
%%
