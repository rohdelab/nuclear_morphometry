clc
clear all
close all

%%
p0=pwd; cd ..; p1=pwd; cd(p0);
sel=16; cl=4;


%%
inp=[p1 '/DATA/data1a/image/']; nm='data1a';
% inp=[p1 '/DATA/data1b/image/']; nm='data2a';
% inp=[p1 '/DATA/data1c/image/']; nm='data3a';
% inp=[p1 '/DATA/data1d/image/']; nm='data4a';

%%
% load([inp 'Img_TOF.mat']);
% cls=unique(label);
% 
% xx1=xx(:,:,find(label==cls(1))); xx1=xx1(:,:,randsample(1:size(xx1,3),16));
% xx2=xx(:,:,find(label==cls(2))); xx2=xx2(:,:,randsample(1:size(xx2,3),16));
% 
% cutt=70;
% 
% im1=[]; cnt=0;
% for a=1:cl
%     im1tmp=[];
%     for b=1:sel/cl
%         cnt=cnt+1;
%         tt=xx1(:,:,cnt); tt=tt(1+cutt:300-cutt,1+cutt:300-cutt);
%         im1tmp=[im1tmp tt];
%     end
%     im1=[im1;im1tmp];
% end
% 
% im2=[]; cnt=0;
% for a=1:cl
%     im2tmp=[];
%     for b=1:sel/cl
%         cnt=cnt+1;
%         tt=xx2(:,:,cnt); tt=tt(1+cutt:300-cutt,1+cutt:300-cutt);
%         im2tmp=[im2tmp tt];
%     end
%     im2=[im2;im2tmp];
% end
% 
% im=[im1 im2];
% [jm] = ib2w(im); 
% 
% save(['Fig_nuclei_' nm '.mat'],'im','jm');

%%
load(['Fig_nuclei_' nm '.mat']);

%%
jm=im; % [jm] = ib2w(im,1); 

jm(481:end,:)=[];
jm(:,1121:end)=[]; jm(:,481:640)=[];

%%
f1=figure; movegui(f1,'northwest')
imshow(im,[]); colormap gray; axis off
pause(1)

f2=figure; movegui(f2,'northeast')
imshow(jm,[]); colormap gray; axis off



