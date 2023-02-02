clc
clear all
close all

warning('off','all')

%%
% tag=1; subtagV={'a','b','c','d'};
tag=3; subtagV={'a'};

%%
p0=pwd; cd ..; p=pwd; cd(p0);

%%
for ii=1:length(subtagV)

    subtag=subtagV{ii};
    disp(['Data - ' num2str(tag) subtag ' started']);

    inp=[p '/DATA/data' num2str(tag) subtag '/'];
    outp=[p '/DATA/' 'P3_CM_WND_dt' num2str(tag) '/DATA/data' num2str(tag) '/org/wndchrm/data' num2str(tag) subtag '/'];

    %%
    load([inp 'image/Img_TOF.mat'])

    if exist(outp)
    else
        mkdir(outp)
    end
    for b=1:size(xx,3)
        nm=[outp 'I' num2str(b) '.tiff'];
        imwrite(uint8(mat2gray(xx(:,:,b))*255),nm);

        disp(['Data - ' num2str(tag) subtag ', image - ' num2str(b)]);
    end
end

%% Next, perform Wnd-charm feature calculations in a local workstation: create wnd_feat.txt