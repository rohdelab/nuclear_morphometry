clc
clear all
close all

%% DATA SELECT
tag=0; subtag=[]; part_type=0; % toy data
tag=0; subtag='z'; part_type=0; % toy data - 2
tag=0; subtag='y'; part_type=0; % toy data - 3

tag=1; subtag='a'; cN=[1277 1689]; part_type=0; % 2966 samples - liver
tag=1; subtag='b'; cN=[161 262]; part_type=0; % 423 samples - thyroid
tag=1; subtag='c'; cN=[590 490]; part_type=0; % 1080 samples - mesothelioma
tag=1; subtag='d'; cN=[5189 6353]; part_type=0; % 11542 samples - melanoma
tag=1; subtag='T'; ust={'a','b','c','d'}; szo=[2966 423 1080 11542]; part_type=0; % All data concatenated, liver, thyroid, melonoma, mesothelioma ...

tag=2; subtag='d'; cN=[2400 2400]; part_type=1; % Martial short: Prostate_epithelial 2 class (? samples): ? (?) & ? (?)

tag=3; subtag='a'; cN=[11310 7190]; part_type=2; % Martial Mesothelioma_sep29, sampled 200/class: cancer & reactive;

tag=4; subtag='a'; cN=[267792 36120]; part_type=2; % Martial Mesothelioma_sep29, full data: cancer & reactive;

%%
p0=pwd; cd ..
inp_labnly=[pwd '/DATA/data' num2str(tag) subtag];
mdpth=[pwd '/DATA/METADATA'];

if subtag=='T'
else
    cd(inp_labnly)
    load(['label' num2str(tag) subtag]);
    load(['patient_label' num2str(tag) subtag]);
end

cd(p0)

Nfold=2;



%% SPLIT & SAVE
vn=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
vpn=[mdpth '/' vn];

%%
if subtag=='T'
    disp('subtag is T')
    for a=1:length(ust)
        tmp_vn=['run5_indsplit_data' num2str(tag) ust{a} '_fold' num2str(Nfold)];
        tmp_vpn=[mdpth '/' tmp_vn]; load(tmp_vpn);
        indB{a}=ind;
    end

    i1=indB{1}{1}; i2=indB{1}{2};
    for a=2:length(ust)
        t1=indB{a}{1}; t2=indB{a}{2};
        t1=t1+sum(szo(1:a-1)); t2=t2+sum(szo(1:a-1));
        i1=[i1(:); t1(:)]; i2=[i2(:); t2(:)];
    end
    ind{1}=i1; ind{2}=i2;
    save(vpn,'ind','-v7.3');
    disp('Completed..')

else
    ISthere=exist([vpn '.mat']);

    if ISthere==0
        disp('Calculating and saving...')
        if part_type==1
            ind=[];
            ind{1}=find(label_patient<25); ind{2}=find(label_patient>24);
        elseif part_type==0
            class=unique(label);
            part_lab=[];
            for a=1:length(class)
                tmp=find(label==class(a));
                tmp2=tmp(:);
                bf=0; TT=[];
                while(1)
                    if length(tmp2)<Nfold
                        tmp22=zeros(Nfold,1);
                        tmp22(1:length(tmp2))=tmp2;
                        tmp2=tmp22;
                        bf=1;
                    end
                    tt=randsample(tmp2,Nfold);
                    TT=[TT tt];
                    if bf==1
                        break;
                    end
                    for c=1:length(tt)
                        tmp2(tmp2==tt(c))=[];
                    end
                end
                part_lab=[part_lab TT];
            end
            part_lab=part_lab';

            for a=1:size(part_lab,2)
                tmp=part_lab(:,a);
                tmp(tmp==0)=[];
                ind{a}=tmp;
            end

        elseif part_type==2

            class=unique(label); ind=[]; ptcc=[];
            for a=1:length(class)
                tmp=find(label==class(a));
                upc=unique(label_patient(tmp));
                ptc = part_arr(upc,Nfold);
                ptcc{a}=ptc;
            end
            for b=1:Nfold
                ii=[];
                for a=1:length(class)
                   ptc=ptcc{a};
                    for c=1:length(ptc{b})
                        t=find(label_patient==ptc{b}(c));
                        if isempty(ii)
                            ii=t(:);
                        else
                            ii=[ii; t(:)];
                        end
                    end

                end
                ind{b}=ii;
            end

        end

        save(vpn,'ind','-v7.3');
        disp('Completed..')
    else
        disp('Result already exists, calculation skipped...')
    end
end