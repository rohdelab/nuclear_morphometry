clc
clear all
close all

warning('off','all')

%%
% tag=1; subtag='a'; tr_te_tag=11; type=1; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
tag=1; subtag='b'; tr_te_tag=11; type=1; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; tr_te_tag=11; type=1; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; tr_te_tag=11; type=1; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
% tag=1; subtag='T'; tr_te_tag=11; type=1; % above four data (and label) concatenated, order: liver, thyroid, meso, mela

% tag=501; subtag='_test1'; tr_te_tag=1; type=1; % 7 types of cancer, Kumar et al.; only testing
% tag=501; subtag='_test2'; tr_te_tag=1; type=1; % slot taken; don't run code for this: see README.txt inside ./DATA/
% tag=502; subtag='_test1'; tr_te_tag=1; type=1; % liver data 3 classes other than 2 above class (3924 samples): FHB (539), FNH (2456), & HCA (929)
% tag=503; subtag='_test1'; tr_te_tag=1; type=1; % thyroid data 4 classes other than 2 above class (2405 samples): FA (756), FC (509), FVPC (300), & NG (840)


% tag=2; subtag='d'; tr_te_tag=11; type=1; % Martial short: Prostate_epithelial 2 class (? samples): ? (?) & ? (?)
% tag=2501; subtag='_test1'; tr_te_tag=1; type=1; % Martial short: same as above


% tag=3; subtag='a'; tr_te_tag=11; type=1; % Martial Mesothelioma_sep29, sampled 200/class: cancer & reactive; both training and testing & only testing
% tag=3501; subtag='_test1'; tr_te_tag=1; type=1; % Martial short: same as above

% tag=4; subtag='a'; tr_te_tag=11; type=2; % Martial Mesothelioma_sep29 - same, full
% tag=4501; subtag='_test1'; tr_te_tag=1; type=2; % Martial full: same as above

%%
reg_str='TOF'; scale_flag=0;
% reg_str='TSOF'; scale_flag=1;

disp(tag); disp(subtag); disp(reg_str);

%%
switch tr_te_tag
    case 11 % both training and testing
        switch tag
            case 1
                switch subtag
                    case 'a'
                        name='liver_normal_fhb_fnh_hca_hcc.mat';
                        nm_patient_label='patient_label_liver.mat';
                        lab_sel=[1 5];
                    case 'b'
                        name='thyroid_nl_fa_fc_fvpc_ng_wifc.mat';
                        nm_patient_label='patient_label_thyroid.mat';
                        lab_sel=[1 6];
                    case 'c'
                        name='mesothelioma_benign_malig.mat';
                        nm_patient_label='patient_label_mesothelioma.mat';
                        lab_sel=[1 2];
                    case 'd'
                        name='melanoma_dn_mm.mat';
                        nm_patient_label='patient_label_melanoma.mat';
                        lab_sel=[1 2];
                    case 'T'
                        nameT{1}='liver_normal_fhb_fnh_hca_hcc.mat';
                        nm_patient_labelT{1}='patient_label_liver.mat';
                        lab_selT{1}=[1 5];
                        nameT{2}='thyroid_nl_fa_fc_fvpc_ng_wifc.mat';
                        nm_patient_labelT{2}='patient_label_thyroid.mat';
                        lab_selT{2}=[1 6];
                        nameT{3}='mesothelioma_benign_malig.mat';
                        nm_patient_labelT{3}='patient_label_mesothelioma.mat';
                        lab_selT{3}=[1 2];
                        nameT{4}='melanoma_dn_mm.mat';
                        nm_patient_labelT{4}='patient_label_melanoma.mat';
                        lab_selT{4}=[1 2];
                    otherwise
                        disp('ERROR!!');
                end
            case 2
                switch subtag
                    case 'a'
                    case 'b'
                    case 'c'
                    case 'd'
                        name='UBC_Prostate_Epithelial_short.mat';
                        nm_patient_label='./UBC_Database/patient_label_Prostate_Epithelial_short.mat';
                        lab_sel=[0 1];
                    otherwise
                        disp('ERROR!!')
                end
            case 3
                switch subtag
                    case 'a'
                        name='UBC_meso_cancer_reactive.mat';
                        nm_patient_label='./UBC_Database/patient_label_meso_cancer_reactive.mat';
                        lab_sel=[0 1];
                    case 'b'
                    case 'c'
                    case 'd'
                    otherwise
                        disp('ERROR!!')
                end
            case 4
                switch subtag
                    case 'a'
                        name='MD/UBC_meso_can_reac_full/';
                        nm_patient_label='./UBC_Database/patient_label_meso_can_reac_full.mat';
                        lab_sel=[0 1];
                    case 'b'
                    case 'c'
                    case 'd'
                    otherwise
                        disp('ERROR!!')
                end
        end
    case 1 % only testing
        switch tag
            case 501
                name='data501_test.mat';
                lab_sel=1:30;
            case 502
                name='liver_normal_fhb_fnh_hca_hcc.mat';
                lab_sel=[2,3,4];
            case 503
                name='thyroid_nl_fa_fc_fvpc_ng_wifc.mat';
                lab_sel=[2,3,4,5];
            case 2501
                name='UBC_Prostate_Epithelial_short.mat';
                lab_sel=[0 1];
            case 3501
                name='UBC_meso_cancer_reactive.mat';
                lab_sel=[0 1];
            case 4501
                name='MD/UBC_meso_can_reac_full/';
                lab_sel=[0 1];
            otherwise
                disp('ERROR!!')
        end
end

%%

p0=pwd;
cd ..
inp=[pwd '/DataBase/DataBase_resized_registered'];
inp_patient_label=[pwd '/DataBase'];
outp=[pwd '/DATA/data' num2str(tag) subtag '/image'];
outp_labnly=[pwd '/DATA/data' num2str(tag) subtag];
mdpth=[pwd '/DATA/METADATA/dta' num2str(tag) subtag '/'];
if exist(mdpth)
else
    mkdir(mdpth)
end


%%
switch tr_te_tag
    case 11 % both training and testing
        if subtag=='T'
            xxT=[]; labelT=[]; label_patientT=[];
            for b=1:length(nameT)
                name=nameT{b}; nm_patient_label=nm_patient_labelT{b}; lab_sel=lab_selT{b};
                cd(inp)
                load([reg_str '_' name]);
                cd(inp_patient_label)
                load(nm_patient_label);
                cd(p0);

                x0=xx; l0=label;
                imgs=xx; labs=label(:); pat_lab=patient_label(:); clear xx; clear label; patient_label;
                xx=[]; label=[]; label_patient=[];
                for a=1:length(lab_sel)
                    xx=cat(3,xx,imgs(:,:,find(lab_sel(a)==labs)));
                    label=[label;labs(find(lab_sel(a)==labs))];
                    label_patient=[label_patient;pat_lab(find(lab_sel(a)==labs))];
                    clc
                    a
                    b
                end

                %
                if scale_flag==0
                    disp('Normalizing (rescaling) images to make the foreground of datasets have same average area. Please wait ...');
                    xx =  avg_area_normalize(xx);
                end
                %

                xxT=cat(3,xxT,xx); labelT=[labelT;label(:)]; label_patientT=[label_patientT;label_patient(:)];
            end
            labelT(labelT>2)=2;
            xx=xxT; label=labelT; label_patient=label_patientT;

            cd(outp)
            save(['Img_' reg_str],'xx','label','-v7.3');

            cd(outp_labnly)
            save(['label' num2str(tag) subtag],'label','-v7.3');
            %     save(['patient_label' num2str(tag) subtag],'label_patient','-v7.3');

            cd(p0)
        else

            if type==1
                cd(inp)
                load([reg_str '_' name]);
                cd(p0);
            elseif type==2
                inp_tmp=[inp '/' name];
            end
            % % %



            if tag ==2
                xx1=xx(:,:,find(label==0)); xx2=xx(:,:,find(label==1));

                xx1_tr=xx1(:,:,1:125); xx2_tr=xx2(:,:,1:125);
                xx1_te=xx1(:,:,126:end); xx2_te=xx2(:,:,126:end);

                xx=[]; label=[]; patient_label=[];
                pcnt=0;

                xx1=xx1_tr; xx2=xx2_tr;
                for rind=1:12
                    pcnt=pcnt+1;
                    if scale_flag==0
                        iind1=randsample(1:size(xx1,3),100);
                        save([mdpth 'iind1_' num2str(pcnt) '.mat'],'iind1')
                    else
                        load([mdpth 'iind1_' num2str(pcnt) '.mat'])
                    end
                    x1=xx1(:,:,iind1); l1=zeros(100,1); p1=pcnt*ones(100,1);
                    xx=cat(3,xx,x1); label=[label;l1]; patient_label=[patient_label;p1];

                    pcnt=pcnt+1;
                    if scale_flag==0
                        iind2=randsample(1:size(xx2,3),100);
                        save([mdpth 'iind2_' num2str(pcnt) '.mat'],'iind2')
                    else
                        load([mdpth 'iind2_' num2str(pcnt) '.mat'])
                    end
                    x2=xx2(:,:,iind2); l2=ones(100,1); p2=pcnt*ones(100,1);
                    xx=cat(3,xx,x2); label=[label;l2]; patient_label=[patient_label;p2];
                end

                xx1=xx1_te; xx2=xx2_te;
                for rind=1:12
                    pcnt=pcnt+1;
                    if scale_flag==0
                        iind1=randsample(1:size(xx1,3),100);
                        save([mdpth 'iind1_' num2str(pcnt) '.mat'],'iind1')
                    else
                        load([mdpth 'iind1_' num2str(pcnt) '.mat'])
                    end
                    x1=xx1(:,:,iind1); l1=zeros(100,1); p1=pcnt*ones(100,1);
                    xx=cat(3,xx,x1); label=[label;l1]; patient_label=[patient_label;p1];

                    pcnt=pcnt+1;
                    if scale_flag==0
                        iind2=randsample(1:size(xx2,3),100);
                        save([mdpth 'iind2_' num2str(pcnt) '.mat'],'iind2')
                    else
                        load([mdpth 'iind2_' num2str(pcnt) '.mat'])
                    end
                    x2=xx2(:,:,iind2); l2=ones(100,1); p2=pcnt*ones(100,1);
                    xx=cat(3,xx,x2); label=[label;l2]; patient_label=[patient_label;p2];
                end

                cell_label=label;
            else

                % only for tag 4?

                cd(inp_patient_label)
                load(nm_patient_label);
                pat_Num=length(unique(patient_label));
                cd(p0);
            end



            % % %
            if type==1
                x0=xx; l0=label;
                imgs=xx; labs=label(:); pat_lab=patient_label(:); clear xx; clear label; patient_label;
                xx=[]; label=[]; label_patient=[];
                for a=1:length(lab_sel)
                    xx=cat(3,xx,imgs(:,:,find(lab_sel(a)==labs)));
                    label=[label;labs(find(lab_sel(a)==labs))];
                    label_patient=[label_patient;pat_lab(find(lab_sel(a)==labs))];
                    clc
                    a
                end

                %
                if scale_flag==0
                    disp('Normalizing (rescaling) images to make the foreground of datasets have same average area. Please wait ...');
                    xx =  avg_area_normalize(xx);
                end
                %

                if exist(outp)
                else
                    mkdir(outp);
                end
                cd(outp)
                save(['Img_' reg_str],'xx','label','-v7.3');

                cd(outp_labnly)
                save(['label' num2str(tag) subtag],'label','-v7.3');
                save(['patient_label' num2str(tag) subtag],'label_patient','-v7.3');

            elseif type==2

                

                cd(inp_tmp); DD=dir('*.mat'); cd(p0);

                label=[]; label_patient=[];
                for patind=1:pat_Num % length(DD)

                    %                     disp(['Patient: ' DD(patind).name]);
                    %                     load([inp_tmp DD(patind).name]);

                    disp(['Patient: ' num2str(patind)]); nm=[reg_str '_pat_dt' num2str(patind) '.mat'];
                    load([inp_tmp nm]);

                    %
                    if scale_flag==0
                        disp('Normalizing (rescaling) images to make the foreground of datasets have same average area. Please wait ...');
                        xx =  avg_area_normalize(xx);
                    end
                    %

                    if exist(outp)
                    else
                        mkdir(outp);
                    end
                    cd(outp)

                    label=[label; ls(:)]; label_patient=[label_patient; pls(:)];

                    lab=ls; p_lab=pls;
                    %                     save(['Img_' DD(patind).name],'xx','lab','p_lab','-v7.3');
                    save(['Img_' nm],'xx','lab','p_lab','-v7.3');



                    cd(p0)
                end

                cd(outp_labnly)
                save(['label' num2str(tag) subtag],'label','-v7.3');
                save(['patient_label' num2str(tag) subtag],'label_patient','-v7.3');
                cd(p0)

            end


        end
    case 1 % only testing


        if type==1
            cd(inp)
            load([reg_str '_' name]);
            cd(p0);
        elseif type==2
            inp_tmp=[inp '/' name];
        end


        if type==1

            x0=xx; l0=label;
            imgs=xx; labs=label(:); clear xx; clear label;


            xx=[]; label=[];
            for a=1:length(lab_sel)
                xx=cat(3,xx,imgs(:,:,find(lab_sel(a)==labs)));
                label=[label;labs(find(lab_sel(a)==labs))];
                clc
                a
            end
            %         xx=imgs; label=labs;

            %
            if scale_flag==0
                disp('Normalizing (rescaling) images to make the foreground of datasets have same average area. Please wait ...');
                xx =  avg_area_normalize(xx);
            end
            %

            if exist(outp)
            else
                mkdir(outp);
            end
            cd(outp)
            save(['Img_' reg_str],'xx','label','-v7.3');

            cd(outp_labnly)
            save(['label' num2str(tag) subtag],'label','-v7.3');

        elseif type==2

            if tag==4501
                    pat_Num=96;
            end

            

            cd(inp_tmp); DD=dir('*.mat'); cd(p0);

            label=[]; label_patient=[];
            for patind=1:pat_Num % length(DD)

                %                 disp(['Patient: ' DD(patind).name]);
                %                 load([inp_tmp DD(patind).name]);

                disp(['Patient: ' num2str(patind)]); nm=[reg_str '_pat_dt' num2str(patind) '.mat'];
                load([inp_tmp nm]);

                %
                if scale_flag==0
                    disp('Normalizing (rescaling) images to make the foreground of datasets have same average area. Please wait ...');
                    xx =  avg_area_normalize(xx);
                end
                %
                if exist(outp)
                else
                    mkdir(outp);
                end
                cd(outp)

                label=[label; ls(:)]; label_patient=[label_patient; pls(:)];

                lab=ls; p_lab=pls;
                %                 save(['Img_' DD(patind).name],'xx','lab','p_lab','-v7.3');
                save(['Img_' nm],'xx','lab','p_lab','-v7.3');

                cd(p0)

            end
            cd(outp_labnly)
            save(['label' num2str(tag) subtag],'label','-v7.3');

            cd(p0)
        end

        cd(p0)
end


