clc
clear all
close all

%%
tag=0; subtag=[]; % toy data

% tag=1; subtag='a'; cN=[1277 1689]; % 2966 samples - liver
% tag=1; subtag='b'; cN=[161 262]; % 423 samples - thyroid
% tag=1; subtag='c'; cN=[590 490]; % 1080 samples - mesothelioma
% tag=1; subtag='d'; cN=[5189 6353]; % 11542 samples - melanoma

%%
p0=pwd; cd ..
mdpth=[pwd '/DATA/METADATA'];
ldp=[pwd '/DATA/data' num2str(tag) subtag '/rcdt'];

cd(ldp)
load Rcdt
xx=u; I0=I0; label=label;

cd(p0)

Nfold=10;

%%
vn=['run5_indsplit_data' num2str(tag) subtag '_fold' num2str(Nfold)];
vpn=[mdpth '/' vn];

load(vpn)

%%
fold_tr=[1 3 4 5 7 8 9 10]; fold_te=[2 6];

[xx_tr,label_tr,xx_te,label_te]=fn_trte_split(xx,label,ind,fold_tr,fold_te);
