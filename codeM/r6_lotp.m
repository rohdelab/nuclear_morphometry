clc
clear all
% close all

warning('off','all')

%%
% tag=0; subtag=[]; % toy data % subtag is a character like 'a', 'b', etc.
% tag=0; subtag='z'; % toy data - 2
% tag=0; subtag='y'; % toy data - 3

% tag=1; subtag='a'; tr_te_tag=11; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; tr_te_tag=11; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; tr_te_tag=11; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; tr_te_tag=11; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
% tag=1; subtag='T'; tr_te_tag=11; % All data concatenated, liver, thyroid, melonoma, mesothelioma ...


% tag=501; subtag='_test1'; tr_te_tag=1; train_tag=1; % 7 types of cancer, Kumar et al.; only testing
% tag=501; subtag='_test2'; tr_te_tag=1; % slot taken; don't run code for this: see README.txt inside ./DATA/
% tag=502; subtag='_test1'; tr_te_tag=1; train_tag=1; % liver data 3 classes other than 2 above class (3924 samples): FHB (539), FNH (2456), & HCA (929)
% tag=503; subtag='_test1'; tr_te_tag=1; train_tag=1; % thyroid data 4 classes other than 2 above class (2405 samples): FA (756), FC (509), FVPC (300), & NG (840)


% tag=2; subtag='d'; tr_te_tag=11; % Martial short: Prostate_epithelial 2 class (? samples): ? (?) & ? (?)
% tag=2501; subtag='_test1'; tr_te_tag=1; train_tag=1; % Martial short: same as above


% tag=3; subtag='a'; tr_te_tag=11; % Martial Mesothelioma_sep29, sampled 200/class: cancer & reactive; both training and testing & only testing
tag=3501; subtag='_test1'; tr_te_tag=1; train_tag=1; % Martial short: same as above


%%
reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'
reg_str_both=1; % 0, 1
Nms=500; %1000; %500;
paral=1;
I0_seed=21;

ref_select=1; % 1 2 3 4. 1 for none; 2, 3 for nothing; 4 for r0

%%
p0=pwd;
cd ..
inp=[pwd '/DATA/data' num2str(tag) subtag '/image'];
outp=[pwd '/DATA/data' num2str(tag) subtag '/lotp'];


switch tr_te_tag
    case 11 % both training and testing
        inpI0=[pwd '/DATA/data' num2str(tag) 'T'];
        cd([inpI0 '/image'])
        if reg_str_both==1
            load(['Img_' 'TSOF'  '.mat']); % TSOF used as reference for both TOF and TSOF
        else
            load(['Img_' reg_str  '.mat']);
        end
        cd(p0);
        imgsall=double(xx); clear xx
    case 1 % only testing
        inpI0=[pwd '/DATA/data' num2str(train_tag) 'T'];
        cd([inpI0 '/image'])
        if reg_str_both==1
            load(['Img_' 'TSOF'  '.mat']); % TSOF used as reference for both TOF and TSOF
        else
            load(['Img_' reg_str  '.mat']);
        end
        cd(p0);
        imgsall=double(xx); clear xx
end
for i=1:size(imgsall,3)
    imgsall(:,:,i)=imgsall(:,:,i)/sum(vec(imgsall(:,:,i)));
end
I0=mean(imgsall,3);

switch ref_select
    case 1
    case 2
        [X,Y]=meshgrid(linspace(-5,5,length(I0)),linspace(-5,5,length(I0)));
        I0=double(X.^2+Y.^2<3);
        I0=I0/sum(vec(I0));
    case 3
        tmp=gaussian1d(linspace(-5,5,length(I0)),1,0); tmp=tmp(:);
        I0=tmp*tmp';
        I0=I0/sum(vec(I0));
    case 4
        [X,Y]=meshgrid(linspace(-5,5,length(I0)),linspace(-5,5,length(I0)));
        msk=double(X.^2+Y.^2<1.5);
        tmp=gaussian1d(linspace(-5,5,length(I0)),1,0); tmp=tmp(:);
        I0=tmp*tmp'; I0=I0.*msk;
        I0=I0/sum(vec(I0));
        ref_str='r0';
    otherwise
end

cd(inp)
load(['Img_' reg_str  '.mat']);
cd(p0);
imgs=double(xx); clear xx
for i=1:size(imgs,3)
    imgs(:,:,i)=imgs(:,:,i)/sum(vec(imgs(:,:,i)));
end


tm=tic;
clc; display(['Calculating LOTP Embedding - ' num2str(tag) subtag ' - ' reg_str]);
[Pl,P]=particleApproximation(imgs,Nms,paral);
rng(I0_seed); [Pl_tem,P_tem]=img2pts_Lloyd(I0,Nms);
[ptcl_wght,LOT_coord,var1]=LOT_LinearEmb(P_tem,Pl_tem,P,Pl,paral);

for a=1:size(LOT_coord,2)
    u(:,a)=reshape((LOT_coord{a})',2*size(ptcl_wght,2),1);
end
toc(tm)

if exist(outp)
else
    mkdir(outp);
end
cd(outp)
switch ref_select
    case 1
        save(['Lotp_' reg_str],'u','ptcl_wght','label','-v7.3');
    case 2
    case 3
    case 4
        save(['Lotp_' ref_str '_' reg_str],'u','ptcl_wght','label','-v7.3');
    otherwise
end

cd(p0)

d=sprintf(['\nFinished r6_lotp_END.m: loaded image space data, applied particle lot (lotp) transform, saved lotp space data - ' num2str(tag) subtag '\n']); disp(d);






















%%                              PAST VERSION -- SAVEV0 version
% % % clc
% % % clear all
% % % close all
% % %
% % % warning('off','all')
% % %
% % % %%
% % % % tag=0; subtag=[]; % toy data % subtag is a character like 'a', 'b', etc.
% % % % tag=0; subtag='z'; % toy data - 2
% % % % tag=0; subtag='y'; % toy data - 3
% % %
% % % % tag=1; subtag='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% % % tag=1; subtag='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% % % % tag=1; subtag='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% % % % tag=1; subtag='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)
% % % % tag=1; subtag='T'; % All data concatenated, liver, thyroid, melonoma, mesothelioma ...
% % %
% % % %%
% % % reg_str='TOF';
% % % Nms=500; %1000; %500;
% % % paral=0;
% % %
% % % %%
% % % p0=pwd;
% % % cd ..
% % % inp=[pwd '/DATA/data' num2str(tag) subtag '/image'];
% % % inpI0=[pwd '/DATA/data' num2str(tag) 'T'];
% % % outp=[pwd '/DATA/data' num2str(tag) subtag '/lotp'];
% % %
% % % cd(inp)
% % % load(['Img_' reg_str  '.mat']);
% % %
% % % cd(p0);
% % % imgs=mat2gray(xx); clear xx;
% % % for i=1:size(imgs,3)
% % %     imgs(:,:,i)=imgs(:,:,i)/sum(vec(imgs(:,:,i)));
% % % end
% % %
% % % if subtag=='T'
% % %     load([inpI0 '/I0_template']) % Using a common template as per Soheil's question
% % % else
% % %     I0=mean(imgs,3);
% % % end
% % %
% % % clc; display(['Calculating LOTP Embedding']);
% % % [Pl,P]=particleApproximation(imgs,Nms,paral);
% % % [Pl_tem,P_tem]=img2pts_Lloyd(I0,Nms);
% % % [ptcl_wght,LOT_coord,var1]=LOT_LinearEmb(P_tem,Pl_tem,P,Pl,paral);
% % %
% % % for a=1:size(LOT_coord,2)
% % %     u(:,a)=reshape((LOT_coord{a})',2*size(ptcl_wght,2),1);
% % % end
% % %
% % % cd(outp)
% % % save(['Lotp_' reg_str],'u','ptcl_wght','label','-v7.3');
% % %
% % % cd(p0)
% % %
% % % d=sprintf('\nFinished r6_lotp_END.m: loaded image space data, applied particle lot (lotp) transform, saved lotp space data.\n'); disp(d);
