clc
clear all
close all

warning('off','all')

%%
% tag=1; subV={'a','b','c','d'};
tag=3; subV={'a'};

%%
p0=pwd; cd ..; p=pwd; cd(p0);

inp=[p '/DATA/WND_receive/data' num2str(tag) '/'];
outp=[p '/DATA/data' num2str(tag)];

%%
inp2=[p '/DATA/data' num2str(tag)];
labV=[]; N=0;
for a=1:length(subV)
    subtag=subV{a};
    inp3=[inp2 subtag '/'];
    load([inp3 'label' num2str(tag) subtag '.mat']);
    labV{a}=label;
    N=N+length(label);
end

%%
f=fopen([inp 'wnd_feat.txt'],'r');

all_data=[];
for a=1:N
    vect=fscanf(f,'%f'); vect(end)=[];
    temp=fscanf(f,'%s/n')
    i=find(temp=='/'); i1=i(end-1); i2=i(end);

    sbt=temp(i2-1); im_no=str2num(temp(i2+2:end-5));

    if tag==1
        switch sbt
            case 'a'
                all_data{1}(:,im_no)=vect(:);
            case 'b'
                all_data{2}(:,im_no)=vect(:);
            case 'c'
                all_data{3}(:,im_no)=vect(:);
            case 'd'
                all_data{4}(:,im_no)=vect(:);
            otherwise
        end
    elseif tag==3
        switch sbt
            case 'a'
                all_data{1}(:,im_no)=vect(:);
            otherwise
        end
    end
    disp([num2str(a) ' of ' num2str(N) ' completed..: ' num2str(a*100/N) ' %']);
end
fclose(f);

%%
for a=1:length(subV)
    subtag=subV{a};
    path1=[outp subtag '/nfe/'];
    if exist(path1)
    else
        mkdir(path1)
    end
    u=all_data{a};
    label=labV{a};
    save([path1 'Nfe_TOF.mat'],'u','label','-v7.3');
end

labV
all_data
