clc
clear all
close all
warning('off','all')

%% PARAMETER SELECT
TagVec = 1; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=5; Nfold=2;
reg_str='TOF'; % 'TOF' 'TSOF'

%% BEGIN
for i_tag=TagVec
    clearvars -except i_tag TagVec Ndir Nfold reg_str
    %% SELECT DATA
    switch i_tag
        case 0
            tag=0; 
            subtag{1}='y'; % toy data - 3
            subtag{2}='z'; % toy data - 2
            subtag{3}=[]; % toy data % subtag is a character like 'a', 'b', etc.
            
            nme=['D_composite_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
        case 1
            tag=1; 
            subtag{1}='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
            subtag{2}='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
            subtag{3}='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
            subtag{4}='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
            
            nme=['D_composite_' num2str(tag)];
            ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'};
        otherwise
            disp('Terminating...')
    end
    %%
    p0=pwd; cd ..; pp=pwd;
    mdpth=[pwd '/DATA/METADATA'];
    respth=[pwd '/RESULTS/EXP05_ex2'];
    Exp_str = 'Exp05_ex2'; Exp_strF = '/EXP05_ex2/';
    cd(p0);
    CL={[0.93 0.45 0.10],[0.09 0.44 0.65],[0.27 0.62 0.25],[0.79 0.15 0.18],[0.4 0.4 0.4]};
    
    %%
    wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
    for foldset=1:Nfold
        vnm1=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
        vnm2=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
        vnm3=[Exp_str '_data' num2str(tag) subtag{1} subtag{2} '_reg' reg_str '___fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
         
        vnm1pnm=[mdpth Exp_strF vnm1];
        vnm2pnm=[mdpth Exp_strF vnm2];
        vnm3pnm=[mdpth Exp_strF vnm3];
        
        ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
        
        if ISthere1==0|ISthere2==0|ISthere3==0
            disp('ERROR!! Result doesnot exist...'); break;
        else
            disp([Exp_str ' testB on Data-' num2str(tag) subtag{1} subtag{2} subtag{3} subtag{4} '__: Calculating and saving... Fold ' num2str(foldset) ' of ' num2str(Nfold)]); disp(' ');
            fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 2-fold split
            
            
            load(vnm2pnm);
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
                
                [~,~,xx_te,lab_te]=fn_trte_split2(data_te,label_te,ind,fold_tr,fold_te);
                [projPLDA_te,axis_image,dm,dt,dt_tr]=proj_n_calc_PLDAl(xx_te,lab_te,PLDA_directions,viz_plda,projPLDA_tr{a},lab_tr{a});
                DMB(:,a,foldset)=dm(:);
                DTB(a,foldset)=dt;
                DT_trB(a,foldset)=dt_tr;
                save(vnm4pnm,'projPLDA_te','-v7.3');
                
            end
            
            f01=figure('position',[0 0 400*size(axis_image,1)/size(axis_image,2) 400]);
            ax1=axes; imagesc(flipud(axis_image'));
            xticks(linspace(1+(length(axis_image)/(2*size(PLDA_directions,2))),length(axis_image)-(length(axis_image)/(2*size(PLDA_directions,2))),size(PLDA_directions,2)));
            xticklabels(1:size(PLDA_directions,2)); yticks([]);
            colormap(ax1,'gray');  xlabel('Modes of variations','interpreter','latex');
            ylabel({'Malignancy','\bf{------------------------$\mathbf{\to}$}'},'interpreter','latex');
            colorbar(ax1,'westoutside','Ticks',linspace(min(axis_image(:)),max(axis_image(:)),5),'TickLabels',{'Much less','Less','-- Neutral --','More','Even more'});
            set(ax1,'fontsize',15);
            title(['Composite model$^{~(' num2str(foldset) '~of~' num2str(Nfold) ')}$'],'interpreter','latex','fontsize',18);
            
            saveas(f01,[respth '/' nme '/testB/axis_set' num2str(foldset) '_' reg_str '.png']);
            
            f02=figure('position',[0 0 1200 450]);
            subplot 121; hold on;
            for a=1:size(DMB,2)
                plot(DMB(:,a,foldset),'o--','markersize',9,'color',CL{a},'markerfacecolor',CL{a});
            end
            grid on;
            plot([0,1+size(DMB,1)],[0 0],'k');
            xticks(1:size(PLDA_directions,2)); xtickangle(45);
            xlim([0 1+size(PLDA_directions,2)]);
            xlabel('Modes of variations','interpreter','latex');
            ylabel('Migration of the mean','interpreter','latex');
            set(gca,'fontsize',15);
            title(['Mean deviation - composite model'],'interpreter','latex','fontsize',18);
            subplot 122; hold on;
            for a=1:size(DMB,2)
                stem((size(DMB,2)-(a-1))*sign(DMB(:,a,foldset)),'o--','markersize',9,'color',CL{a},'markerfacecolor',CL{a});
            end
            plot([0,1+size(DMB,1)],[0 0],'k'); grid on;
            text(0.65,size(DMB,2)+0.7,'Correct direction','fontsize',12,'color',[.4 .4 .4]);
            text(0.65,-(size(DMB,2)+0.7),'Wrong direction','fontsize',12,'color',[.4 .4 .4]);
            xticks(1:size(PLDA_directions,2)); xtickangle(45);
            xlim([0 1+size(PLDA_directions,2)]);
            ylim([-(size(DMB,2)+1) (size(DMB,2)+1)]);
            xlabel('Modes of variations','interpreter','latex');
            ylabel('Migration of the mean (quantized)','interpreter','latex');
            legend([otnme],'location','southoutside');
            set(gca,'fontsize',15);
            title(['Model - ' num2str(foldset) ' of ' num2str(Nfold)],'interpreter','latex','fontsize',18)
            
            saveas(f02,[respth '/' nme '/testB/meanMig_set' num2str(foldset) '_' reg_str '.png']);
            
            f03=figure('position',[0 0 600 350]);
            bar([DTB(:,foldset) DT_trB(:,foldset)]); grid on;
            legend('Between class','Within class');
            xticklabels(otnme); set(gca,'fontsize',15);
            title(['Wasserstein distance - composite model'],'interpreter','latex','fontsize',18);
            
            saveas(f03,[respth '/' nme '/testB/wassDist_set' num2str(foldset) '_' reg_str '.png']);
            
            pause(10); close all;
        end
        waitbar(foldset/Nfold,wh);
    end
    close(wh); close all;
    
    n1=size(DMB,1); n2=size(DMB,2); n3=size(DMB,3);
    D=[];
    for a=1:n1
        for b=1:n2
            for c=1:n3
                D(c,b,a)=DMB(a,b,c);
            end
        end
    end
    DD=[];
    for a=1:size(D,3)
        temp=D(:,:,a);
        DD=[DD temp];
        DD(:,size(DD,2)+1)=nan;
    end
    str=['MeanMigAll_' nme(3:5)];
    for a=1:length(otnme)
        str=[str '_' otnme{a}(1:3)];
    end
    str=[str '_' reg_str '.csv'];
    csvwrite([respth '/' nme '/testB/' str],DD)
    d=sprintf(['\nFinished ' Exp_str '-testB on Data' num2str(tag) subtag{1} subtag{2} '__: projected test data on PLDA directions in LOTP space; saved PLDA projections; saved the plots of the migration of the mean and axis visualization figures; saved the numerical results in excel.\n\n']); disp(d);
end
%%
