clc
clear all
close all

%%
% im=imread(['allseg.png']);
% im=mean(im,3);
% [jm] = ib2w(im,1); 

%%
km=[]; km1=[];
load(['Fig_nuclei_' 'data1a' '.mat']); [jm] = ib2w(im,1); km1=jm;
load(['Fig_nuclei_' 'data2a' '.mat']); [jm] = ib2w(im,1); km1=[km1 jm];
load(['Fig_nuclei_' 'data3a' '.mat']); [jm] = ib2w(im,1); km2=jm;
load(['Fig_nuclei_' 'data4a' '.mat']); [jm] = ib2w(im,1); km2=[km2 jm];
km=[km1;km2];



%%
f1=figure; movegui(f1,'northwest')
imshow(im,[]); colormap gray; axis off
pause(1)

f2=figure; movegui(f2,'northeast')
imshow(jm,[]); colormap gray; axis off

f3=figure; movegui(f3,'north')
imshow(km,[]); colormap gray; axis off



