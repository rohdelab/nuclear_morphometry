clc
clear all
close all

warning('off','all')

%%
% tag=1; subtagV={'a','b','c','d'};
tag=3; subtagV={'a'};

select_type=2; % 1,2,3, ...

%%
p0=pwd; cd ..; p=pwd; cd(p0);

%%
for ii=1:length(subtagV)

    subtag=subtagV{ii};
    disp(['Data - ' num2str(tag) subtag ', NFE type - ' num2str(select_type) ' started']);

    inp=[p '/DATA/data' num2str(tag) subtag '/'];
    outp=[p '/DATA/data' num2str(tag) subtag '/nfe/'];

    %%
    switch select_type
        case 1
            % Type 1
            load([inp 'lotp/Lotp_TOF.mat'])
            u=10*rand(size(u));

        case 2
            % Type 2
            load([inp 'image/Img_TOF.mat'])

            u=[];
            for a=1:size(xx,3) % can do parfor if you want
                [tf] = morh_feat_fn(xx(:,:,a));
                u(:,a)=tf(:);
                disp(['Data - ' num2str(tag) subtag ', NFE type - ' num2str(select_type) '; image - ' num2str(a)])
            end

        otherwise
            % Error
            disp('error !!')
    end

    disp(size(u))
    disp(size(label))

    %%
    if exist(outp)
    else
        mkdir(outp)
    end
    save([outp 'Nfe_TOF_MATLAB.mat'],'u','label','-v7.3');
end
