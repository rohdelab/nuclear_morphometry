clc
clear all
close all

warning('off','all')

%%
tag=0; subtag=[]; nme='D_toy'; ote={'1a','1b','1c','1d'}; otnme={'liver','thyroid','mesothelioma','melanoma'}; % toy data % subtag is a character like 'a', 'b', etc.

% tag=1; subtag='a'; nme='D_liver'; ote={'1b','1c','1d'}; otnme={'thyroid','mesothelioma','melanoma'}; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; nme='D_thyroid'; ote={'1a','1c','1d'}; otnme={'liver','mesothelioma','melanoma'}; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; nme='D_mesothelioma'; ote={'1a','1b','1d'}; otnme={'liver','thyroid','melanoma'}; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; nme='D_melanoma'; ote={'1a','1b','1c'}; otnme={'liver','thyroid','mesothelioma'}; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

%% PARAMETER SELECT
Ndir=20; Nfold=10;
reg_str='TOF';

%%
p0=pwd; cd ..; pp=pwd;
mdpth=[pwd '/DATA/METADATA'];
respth=[pwd '/RESULTS/EXP01'];
Exp_str = 'Exp01'; Exp_strF = '/EXP01/';
cd(p0);

%%
wh=waitbar(0,['Saving ' num2str(Nfold) ' folds...']);
for foldset=1:Nfold
    vnm1=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_PLDA_directions'];
    vnm2=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_tr'];
    vnm3=[Exp_str '_data' num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_viz_plda'];
    
    
    vnm1pnm=[mdpth Exp_strF vnm1];
    vnm2pnm=[mdpth Exp_strF vnm2];
    vnm3pnm=[mdpth Exp_strF vnm3];
    
    
    
    ISthere1=exist([vnm1pnm '.mat']); ISthere2=exist([vnm2pnm '.mat']); ISthere3=exist([vnm3pnm '.mat']);
    
    
    if ISthere1==0|ISthere2==0|ISthere3==0
        disp('ERROR!! Result doesnot exist...'); break;
    else
        disp([Exp_str ' testB on Data ' num2str(tag) subtag ': Calculating and saving... ' num2str(foldset) ' of ' num2str(Nfold)]);
        fold_tr=[1:foldset-1 foldset+1:Nfold]; fold_te=[foldset]; % 10-fold split
        
        vnm4=[Exp_str '_data' num2str(tag) subtag num2str(tag) subtag '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
        vnm4pnm=[mdpth Exp_strF vnm4];
        
        indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
        indpnm=[mdpth '/' indnm];
        
        load(indpnm)
        load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
        
        inp_te=[pp '/DATA/data' num2str(tag) subtag '/rcdt'];
        load([inp_te '/Rcdt_' reg_str]);
        data_te=u; label_te=label;
        
        [~,~,xx_te,lab_te]=fn_trte_split(data_te,label_te,ind,fold_tr,fold_te);
        [projPLDA_te,axis_image,dm]=proj_n_calc_PLDAr(xx_te,lab_te,PLDA_directions,viz_plda);
        DMB(:,1,foldset)=dm(:);
        save(vnm4pnm,'projPLDA_te','-v7.3');
        
        
        for a=1:length(ote)
            vnm4=[Exp_str '_data' num2str(tag) subtag ote{a} '_reg' reg_str '_fold' num2str(foldset) 'of' num2str(Nfold) '_retdir' num2str(Ndir) '_projPLDA_te'];
            vnm4pnm=[mdpth Exp_strF vnm4];
            
            indnm=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
            indpnm=[mdpth '/' indnm];
            
            load(indpnm)
            load(vnm1pnm); load(vnm2pnm); load(vnm3pnm);
            
            inp_te=[pp '/DATA/data' ote{a} '/rcdt'];
            load([inp_te '/Rcdt']);
            data_te=u; label_te=label;
            xx_te=data_te;  lab_te=label_te;
            [projPLDA_te,axis_image,dm]=proj_n_calc_PLDAr(xx_te,lab_te,PLDA_directions,viz_plda);
            DMB(:,1+a,foldset)=dm(:);
            save(vnm4pnm,'projPLDA_te','-v7.3');
            
        end
        
        f01=figure;
        imshow(axis_image',[]);
        xlabel('Axis number'); ylabel('Variation along axis');
        title([nme(3:end) ' model'])
        set(gca,'fontsize',12);
        saveas(f01,[respth '/' nme '/testB/axis_set' num2str(foldset) '_' reg_str '.png']);
        
        f02=figure('position',[0 0 1200 450]);
        subplot 121
        plot(DMB(:,1,foldset),'s--','markersize',9); hold on;
        plot(DMB(:,2:end,foldset),'o--','markersize',9); grid on;
        plot([1,size(DMB,1)],[0 0],'k')
        xlabel('Axis number'); ylabel('Mean migration');
        title([nme(3:end) ' model'])
        legend([nme(3:end),otnme]); set(gca,'fontsize',12);
        subplot 122
        plot(sign(DMB(:,1,foldset)),'s--','markersize',9); hold on;
        for a=2:size(DMB,2)
            plot(a*sign(DMB(:,a,foldset)),'o--','markersize',9);
        end
        plot([1,size(DMB,1)],[0 0],'k'); grid on; 
        xlabel('Axis number'); ylabel('Mean migration');
        title([nme(3:end) ' model'])
        legend([nme(3:end),otnme]); set(gca,'fontsize',12);
        saveas(f02,[respth '/' nme '/testB/meanMig_set' num2str(foldset) '_' reg_str '.png']);
        
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

%%
d=sprintf(['\nFinished ' Exp_str '-testB on Data' num2str(tag) subtag ':??.\n']); disp(d);
