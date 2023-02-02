clc
clear all
% close all

warning('off','all')

%%
% tag=4; subtag='a'; tr_te_tag=11; patNum=96; % Martial Mesothelioma_sep29, full data: cancer & reactive; both training and testing & only testing
tag=4501; subtag='_test1'; tr_te_tag=1; patNum=96; train_tag=1; % Martial full: same as above

%%
reg_str='TOF'; % 'TOF' 'TSOF'
reg_str_both=1; % 0, 1
Nms=500; %1000; %500;
paral=1;
I0_seed=21;

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


%%
donecnt=0;
for ii=1:patNum
    if exist(outp)
    else
        mkdir(outp);
    end
    namei=[inp '/Img_' reg_str '_pat_dt' num2str(ii) '.mat'];
    nameo=[outp '/Lotp_' reg_str '_pat_dt' num2str(ii) '.mat'];


    if exist(nameo)
        disp(['Patient ' num2str(ii) ' results exists.']);
        donecnt=donecnt+1;
    else
        load(namei); label=lab;
        imgs=double(xx); clear xx
        for i=1:size(imgs,3)
            imgs(:,:,i)=imgs(:,:,i)/sum(vec(imgs(:,:,i)));
        end


        clc; display(['Calculating LOTP Embedding - ' num2str(tag) subtag ' - ' reg_str]);
        [Pl,P]=particleApproximation(imgs,Nms,paral);
        rng(I0_seed); [Pl_tem,P_tem]=img2pts_Lloyd(I0,Nms);
        [ptcl_wght,LOT_coord,var1]=LOT_LinearEmb(P_tem,Pl_tem,P,Pl,paral);

        u=[];
        for a=1:size(LOT_coord,2)
            u(:,a)=reshape((LOT_coord{a})',2*size(ptcl_wght,2),1);
        end

        save(nameo,'u','ptcl_wght','label','-v7.3');
    end

end

if donecnt==patNum
    nameo=[outp '/Lotp_' reg_str '.mat'];
    uB=[];lB=[];
    for ii=1:patNum
        namei=[outp '/Lotp_' reg_str '_pat_dt' num2str(ii) '.mat'];
        load(namei);
        uB=[uB u(:,1:length(label(:)))]; lB=[lB;label(:)];
        clc;ii
    end
    u=uB; label=lB;
    save(nameo,'u','ptcl_wght','label','-v7.3');
end












