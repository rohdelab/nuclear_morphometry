clc
clear all
close all

%%
p0=pwd; cd ..
p=pwd;
c=[pwd '/DataBase/DataBase_resized'];
d=[pwd '/DataBase/DataBase_resized_registered'];
cd('./codeM')

%%
name{1}='/liver_normal_fhb_fnh_hca_hcc.mat';
name{2}='/thyroid_nl_fa_fc_fvpc_ng_wifc.mat';
name{3}='/mesothelioma_benign_malig.mat';
name{4}='/melanoma_dn_mm.mat';
tr_te_tag=11; % both training and testing
type=1;

% name='/data501_test.mat';
% tr_te_tag=1; % 7 types of cancer, Kumar et al.; only testing
% type=1;

% name='/UBC_Prostate_Epithelial_short.mat';
% tr_te_tag=111; % Martial Mesothelioma_sep29: cancer & reactive; both training and testing & only testing
% type=1;

% name='/UBC_meso_cancer_reactive.mat';
% tr_te_tag=111; % Martial Mesothelioma_sep29, sampled 200/class: cancer & reactive; both training and testing & only testing
% type=1;

% name='/MD/UBC_meso_can_reac_full/';
% tr_te_tag=111; % Martial Mesothelioma_sep29 - same; take full dataset, not a small portion like 200/class
% type=2;

%%
% Reg_Mode=10; % TO
% Reg_Mode=11; reg_str='TOF'; % TOF
Reg_Mode=15; reg_str='TSOF'; % TSOF

%%
xx=[];

switch tr_te_tag
    case 11 % both training and testing
        for a=1:length(name)
            clc; disp('loading....'); disp(a); disp([name{a} ' - ' reg_str])
            load([c name{a}]); xx=xx_resized; label=label; clear xx_resized
            clc; disp('Converting to grayscale...')
            for b=1:length(label)
                xx(:,:,b)=mat2gray(xx(:,:,b));
            end
            clc; disp([num2str(a) ' - Translation, Scaling, Orientation, Flipping...'])
            xx_r=imageInitialization(xx,Reg_Mode);
            xx=xx_r; clear xx_r
            nm=name{a}; nm2=[nm(1) reg_str '_' nm(2:end)];
            save([d nm2],'xx','label','-v7.3');
        end
    case 1 % only testing
        clc; disp('loading....'); disp([name ' - ' reg_str])
        load([c name]); xx=xx_resized; label=label; clear xx_resized
        clc; disp('Converting to grayscale...')
        for b=1:length(label)
            xx(:,:,b)=mat2gray(xx(:,:,b));
        end
        clc; disp([' - Translation, Scaling, Orientation, Flipping...'])
        xx_r=imageInitialization(xx,Reg_Mode);
        xx=xx_r; clear xx_r
        nm=name; nm2=[nm(1) reg_str '_' nm(2:end)];
        save([d nm2],'xx','label','-v7.3');
    case 111 % both training and testing & only testing
        if type==1
            clc; disp('loading....'); disp([name ' - ' reg_str])
            load([c name]); xx=xx_resized; label=label; clear xx_resized
            clc; disp('Converting to grayscale...')
            for b=1:length(label)
                xx(:,:,b)=mat2gray(xx(:,:,b));
            end
            clc; disp([' - Translation, Scaling, Orientation, Flipping...'])
            xx_r=imageInitialization(xx,Reg_Mode);
            xx=xx_r; clear xx_r
            nm=name; nm2=[nm(1) reg_str '_' nm(2:end)];
            save([d nm2],'xx','label','-v7.3');
        elseif type==2
            clc; disp('loading....'); disp([name ' - ' reg_str])
            ip=[c name];
            cd(ip)
            DD=dir('*.mat');

            cd(p0);

            for b=1:length(DD)
                op=[d name];
                if exist(op)
                else
                    mkdir(op)
                end
                load([ip DD(b).name]); pat=pls(1);

                if exist([op reg_str '_pat_dt' num2str(pat) '.mat'])
                    disp(['Patient - ' num2str(pat) ' - ' reg_str ' exists. Skipping calculation..'])
                else
                    xx=xx_resized; clear xx_resized
                    N=size(xx,3);
                    clc; disp(['Converting to grayscale... - ' num2str(b)])
                    for a=1:N
                        xx(:,:,a)=mat2gray(xx(:,:,a));
                    end
                    disp([' - Translation, Scaling, Orientation, Flipping... - ' num2str(b)])
                    xx_r=imageInitialization(xx,Reg_Mode);
                    xx=xx_r; clear xx_r


                    save([op reg_str '_pat_dt' num2str(pat) '.mat'],'xx','ls','pls','-v7.3');
                end
            end

        end
end
%%