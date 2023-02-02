clc
clear all
close all
% normally r2 then r3, here r3 then r2;
% can follow r2 then r3, or can just run this (r32) and skip both r2 and r3

%%
tag=1;

%%
% Reg_Mode=11; reg_str='TOF_post'; scale_flag=0; % TOF
Reg_Mode=15; reg_str='TSOF_post'; scale_flag=1; % TSOF

%%
switch tag
    case 1
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
        sbtgT={'a','b','c','d'};
    otherwise
        disp('ERROR!!');
end

%%
p0=pwd; cd ..

inp=[pwd '/DataBase/DataBase_resized'];
inp_patient_label=[pwd '/DataBase'];

outp0=[pwd '/DATA/data' num2str(tag)];
outp_labnly0=[pwd '/DATA/data' num2str(tag)];


%%
xxT=[]; labelT=[]; label_patientT=[]; Ncl=[];
for b=1:length(nameT)
    name=nameT{b}; nm_patient_label=nm_patient_labelT{b}; lab_sel=lab_selT{b};
    cd(inp)
    load([name]); xx=xx_resized; label=label; clear xx_resized;
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
    Ncl(b)=length(label);
    
    xxT=cat(3,xxT,xx); labelT=[labelT;label(:)]; label_patientT=[label_patientT;label_patient(:)];
end


clc; disp('Converting to grayscale...')
for a=1:length(labelT)
    xxT(:,:,b)=mat2gray(xxT(:,:,b)); % xxT(:,:,b)/sum(sum(xxT(:,:,b)));
end
clc; disp([num2str(a) ' - Translation, Scaling, Orientation, Flipping...'])
xx_r=imageInitialization(xxT,Reg_Mode);
xxT=xx_r;


Ncl=[0;Ncl(:)];
for a=2:length(Ncl)
    outp=[outp0 sbtgT{a-1} '/image'];
    outp_labnly=[outp_labnly0 sbtgT{a-1}];
    ind=sum(Ncl(1:a-1))+1:sum(Ncl(1:a));
    xx=xxT(:,:,ind); label=labelT(ind); label_patient=label_patientT(ind);
    
    %
    if scale_flag==0
        disp('Normalizing (rescaling) images to make the foreground of datasets have same average area. Please wait ...');
        xx =  avg_area_normalize(xx);
        xxT(:,:,ind)=xx;
    end
    %
    
    cd(outp)
    save(['Img_' reg_str],'xx','label','-v7.3');
    cd(outp_labnly)
    save(['label' num2str(tag) sbtgT{a-1}],'label','-v7.3');
    save(['patient_label' num2str(tag) sbtgT{a-1}],'label_patient','-v7.3');
    cd(p0)
end

outp=[outp0 'T' '/image'];
outp_labnly=[outp_labnly0 'T'];
labelT(labelT>2)=2;
xx=xxT; label=labelT; label_patient=label_patientT;

cd(outp)
save(['Img_' reg_str],'xx','label','-v7.3');
cd(outp_labnly)
save(['label' num2str(tag) 'T'],'label','-v7.3');
cd(p0)


