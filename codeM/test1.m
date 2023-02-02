clc
clear all
close all

%%
train_tag=1;

%%
reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'
Nms=500; %1000; %500;
paral=0;
I0_seed=21;

%%
p0=pwd;
cd ..
% inp=[pwd '/DATA/data' num2str(tag) subtag '/image'];
% outp=[pwd '/DATA/data' num2str(tag) subtag '/lotp'];


inpI0=[pwd '/DATA/data' num2str(train_tag) 'T'];
cd([inpI0 '/image'])
load(['Img_' reg_str  '.mat']);
cd(p0);
imgsall=double(xx); clear xx

for i=1:size(imgsall,3)
    imgsall(:,:,i)=imgsall(:,:,i)/sum(vec(imgsall(:,:,i)));
end
I0=mean(imgsall,3);


clear xx
l=size(I0,1); r1=80; r2=100;
p=ceil(-l/2:l/2-1); q=ceil(-l/2:l/2-1); [P,Q]=meshgrid(p,q);
xx(:,:,1)=P.^2+Q.^2<r1^2; %xx(:,:,2)=P.^2+Q.^2<r2^2;
imgs=double(xx); clear xx



for i=1:size(imgs,3)
    imgs(:,:,i)=imgs(:,:,i)/sum(vec(imgs(:,:,i)));
end


tm=tic;
clc;
[Pl,P]=particleApproximation(imgs,Nms,paral);
rng(I0_seed); [Pl_tem,P_tem]=img2pts_Lloyd(I0,Nms);

% % %
tmp=Pl{1}; fact=0.75;
mb1=mean(tmp(1,:)); mb2=mean(tmp(2,:));
tmp(1,:)=tmp(1,:)-mb1; tmp(2,:)=tmp(2,:)-mb2;
tmp=tmp*fact;
tmp(1,:)=tmp(1,:)+mb1; tmp(2,:)=tmp(2,:)+mb2;
Pl{2}=tmp; P{2}=P{1};

% tmp=Pl_tem; fact1=0.75; fact2=0.5;
% tmp1= tmp; tmp2=tmp;
% mb1=mean(tmp(1,:)); mb2=mean(tmp(2,:));
% 
% tmp1(1,:)=tmp1(1,:)-mb1; tmp1(2,:)=tmp1(2,:)-mb2;
% tmp1=tmp1*fact1;
% tmp1(1,:)=tmp1(1,:)+mb1; tmp1(2,:)=tmp1(2,:)+mb2;
% Pl{1}=tmp1; P{1}=P_tem;
% 
% tmp2(1,:)=tmp2(1,:)-mb1; tmp2(2,:)=tmp2(2,:)-mb2;
% tmp2=tmp2*fact2;
% tmp2(1,:)=tmp2(1,:)+mb1; tmp2(2,:)=tmp2(2,:)+mb2;
% Pl{2}=tmp2; P{2}=P_tem;
% % %


[ptcl_wght,LOT_coord,var1]=LOT_LinearEmb(P_tem,Pl_tem,P,Pl,paral);

for a=1:size(LOT_coord,2)
    u(:,a)=reshape((LOT_coord{a})',2*size(ptcl_wght,2),1);
end
toc(tm)

out=u(:,1)-u(:,2);



