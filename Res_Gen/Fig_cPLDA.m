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

            nme=['Dataset_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
        case 1
            tag=1;
            subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

            nme=['Dataset_' num2str(tag)];

            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};

        otherwise
            disp('Terminating...')
    end
    %%
    p0=pwd; cd ..; pp=pwd;
    addpath([pwd '/codeM'])
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/MAIN'];
    Exp_str = 'Main'; Exp_strF = '/MAIN/';
    cd(p0);
    CL={[0.93 0.45 0.10],[0.09 0.44 0.65],[0.27 0.62 0.25],[0.79 0.15 0.18],[0.4 0.4 0.4]};

    %%
    for foldset=1:1 % Nfold
        if reg_str_both==1
            vnm1=[Exp_str '_data' num2str(tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=[Exp_str '_data' num2str(tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=[Exp_str '_data' num2str(tag) 'ab' '_reg' 'BOTH' '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        else
            vnm1=[Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
            vnm2=[Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
            vnm3=[Exp_str '_data' num2str(tag) 'ab' '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
        end

        vnm1pnm=[mdpth Exp_strF vnm1]; vnm2pnm=[mdpth Exp_strF vnm2]; vnm3pnm=[mdpth Exp_strF vnm3];

        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);

        if ISthere1==0|ISthere2==0|ISthere3==0
            disp('ERROR!! Result doesnot exist...'); break;
        else
            disp([Exp_str ' testB on Data-' num2str(tag) '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split


            load(vnm2pnm);
            SD_spread=3; mvim=7; % set inside proj_n_calc_PLDAl.m - see if matches
            for a=1:length(ote)
                vnm4=[Exp_str '_data' num2str(tag) 'ab' ote{a} '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
                vnm4pnm=[mdpth Exp_strF vnm4];

                indnm=['run5_indsplit_data' ote{a} '_fold' num2str(Nfold)];
                indpnm=[mdpth '/' indnm];
                load(indpnm)

                load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);

                inp_te=[pp '/DATA/data' ote{a} '/lotp'];
                load([inp_te '/Lotp_' reg_str]);
                data_te=u; label_te=label;

                [~,~,xx_te,lab_te]=fn_trte_split2(data_te,label_te,ind,fold_tr,fold_te);
                [projPLDA_te,axis_image,dm,dt,dt_tr]=proj_n_calc_PLDAl(xx_te,lab_te,PLDA_directions,viz_plda,projPLDA_tr{a},lab_tr{a});
                DMB(1:length(dm(:)),a,foldset)=dm(:);
                DTB(a,foldset)=dt;
                DT_trB(a,foldset)=dt_tr;

            end

            s=size(axis_image,2)/mvim; st=s/2; ytc0=st:s:size(axis_image,2); ytc=linspace(-SD_spread,SD_spread,mvim);

            f01=figure('position',[0 0 400 300*size(axis_image,1)/size(axis_image,2)]*1.2);
            subplot 335
            ax1=axes;
            %                 imagesc((axis_image));

            % %
%             [axis_image] = ib2w(axis_image,0);
            % %


            imshow(axis_image); axis on;
            yticks(linspace(1+(length(axis_image)/(2*size(PLDA_directions,2))),length(axis_image)-(length(axis_image)/(2*size(PLDA_directions,2))),size(PLDA_directions,2)));
            yticklabels(1:size(PLDA_directions,2));

            xtc1=[];
            for a=1:length(ytc)
                if a==ceil(length(ytc)/2)
                    xtc1{a}=[num2str(ytc(a))];
                else
                    xtc1{a}=[num2str(ytc(a)) '\sigma'];
                end
            end

            xticks(ytc0); xticklabels(xtc1); xlabel('Projection score'); % yticks([]);
            colormap(ax1,'gray');  ylabel('Modes of variations','interpreter','latex');

            set(ax1,'fontsize',15);
            %             title(['Composite model'],'interpreter','latex','fontsize',18);
            set(gcf,'position',[0 0 400 550]); movegui(f01,'north')

            Fgflnm=['Fig_' 'cPLDA' '/'];
            if exist(Fgflnm)
            else
                mkdir(Fgflnm)
            end
            saveas(f01,[Fgflnm 'composite' '.svg']);
            saveas(f01,[Fgflnm 'composite' '.png']);

            pause(1); close all;
        end
    end

end
%%
