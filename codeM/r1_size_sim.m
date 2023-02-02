clc
clear all
close all

%%
% tag=101; tr_te_tag=11; % Liver all; both training and testing
tag=102; tr_te_tag=11; % Thyroid all; both training and testing
% tag=103; tr_te_tag=11; % Mesothelioma all; both training and testing
% tag=104; tr_te_tag=11; % Melanoma all; both training and testing

% tag=201; tr_te_tag=111; % Martial short: Prostate; both training and testing & only testing

% tag=301; tr_te_tag=111; % Martial Mesothelioma_sep29, sampled 200/class: cancer & reactive; both training and testing & only testing

% tag=401; tr_te_tag=111; % Martial Mesothelioma_sep29 - same; take full dataset, not a small portion like 200/class

% tag=501; tr_te_tag=1; % 7 types of cancer, Kumar et al.; only testing


%%
p0=pwd; cd ..
p=pwd;
c=[pwd '/DataBase'];
c1=[pwd '/Testing_Database/Data' num2str(tag)];
d=[pwd '/DataBase/DataBase_resized'];
cd('./codeM')

new_sz=300; % max image size across all datasets

%%
switch tr_te_tag
    case 11 % both training and testing
        switch tag
            case 101
                path_in=c;
                name='/liver_normal_fhb_fnh_hca_hcc.mat';
                load([c name]);
                label=liver_label;
                xx=liver;
                type=1;
            case 102
                name='/thyroid_nl_fa_fc_fvpc_ng_wifc.mat';
                load([c name]);
                label=thyroid_label;
                xx=thyroid;
                type=1;
            case 103
                name='/mesothelioma_benign_malig.mat';
                load([c name]);
                label=mesothelioma_label;
                xx=mesothelioma;
                type=1;
            case 104
                name='/melanoma_dn_mm.mat';
                load([c name]);
                label=melanoma_label;
                xx=melanoma;
                type=1;
            otherwise
                disp('no');
        end
    case 1 % only testing
        switch tag
            case 501
                path_in=c1;
                name='/data501_test.mat';
                load([c1 name]);
                label=label;
                xx=Img;
                type=1;
            otherwise
                disp('no')
        end
    case 111 % both training and testing & only testing
        switch tag
            case 201
                path_in=[c '/UBC_Database'];
                name='/Prostate_Epithelial_short.mat';
                load([path_in name]);
                label=labs;
                xx=imgs;
                name='/UBC_Prostate_Epithelial_short.mat';
                type=1;
            case 301
                path_in=[c '/UBC_Database'];
                name='/meso_cancer_reactive.mat';
                load([path_in name]);
                label=meso_label;
                xx=meso;
                name='/UBC_meso_cancer_reactive.mat';
                type=1;
            case 401
                path_in=[c '/UBC_Database/MD/meso_can_reac_full/'];
                %                 name='/meso_can_reac_full.mat';
                %                 load([path_in name]);
                %                 label=meso_label;
                %                 xx=meso;
                name='/UBC_meso_can_reac_full';
                type=2;
            otherwise
                disp('no')
        end
    otherwise
        disp('never')
end

%%
if type==1
    clearvars -except xx new_sz d name label
    N=size(xx,3); past_sz=size(xx,1); xx_resized=[];
    for a=1:N
        temp=interp1(linspace(0,1,past_sz),xx(:,:,a),linspace(0,1,new_sz));

        xx_resized(:,:,a)=interp1(linspace(0,1,past_sz),temp',linspace(0,1,new_sz))';

        clc; disp(['Completed: ' num2str(a*100/N) ' %'])
    end

    save([d name],'xx_resized','label','-v7.3');
elseif type==2
    cd(path_in)
    DD=dir('*.mat');
    cd(p0);
    for b=1:length(DD)
        load([path_in DD(b).name]); pat=pls(1);

        N=size(ims,3); past_sz=size(ims,1); xx_resized=[];
        for a=1:N
            temp=interp1(linspace(0,1,past_sz),ims(:,:,a),linspace(0,1,new_sz));

            xx_resized(:,:,a)=interp1(linspace(0,1,past_sz),temp',linspace(0,1,new_sz))';

            clc; disp(['Completed ' num2str(b) ' - ' num2str(a*100/N) ' %'])
        end

        op=[d '/MD' name '/'];
        if exist(op)
        else
            mkdir(op)
        end
        save([op 'pat_dt' num2str(pat) '.mat'],'xx_resized','ls','pls','-v7.3');
    end
end



