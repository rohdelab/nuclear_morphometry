rrclc
clear all
close all

warning('off','all')

%%
% tag=0; subtag=[]; % toy data % subtag is a character like 'a', 'b', etc.
% tag=0; subtag='z'; % toy data - 2
% tag=0; subtag='y'; % toy data - 3

tag=1; subtag='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

%%
reg_str='TOF'; 

%%
p0=pwd;
cd ..
inp=[pwd '/DATA/data' num2str(tag) subtag '/image'];
outp=[pwd '/DATA/data' num2str(tag) subtag '/rcdt'];

cd(inp)
load(['Img_' reg_str '.mat']);

cd(p0);
imgs=mat2gray(xx); clear xx;
for i=1:size(imgs,3)     
    imgs(:,:,i)=imgs(:,:,i)/sum(vec(imgs(:,:,i)));
end

I0=mean(imgs,3);
parfor a=1:size(imgs,3)
    clc; display(['Calculating RCDT Embedding for image # ',num2str(a)])
    [~,~,u(:,:,a),~,~]=RCDT(imgs(:,:,a),I0);
end


cd(outp)
save(['Rcdt_' reg_str],'u','I0','label','-v7.3');

cd(p0)

d=sprintf('\nFinished r5_rcdt.m: loaded image space data, applied rcdt transform, saved rcdt space data.\n'); disp(d);