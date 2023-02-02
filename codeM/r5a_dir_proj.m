clc
clear all
close all

warning('off','all')

%%
tag=0; subtag=[]; % toy data % subtag is a character like 'a', 'b', etc.

% tag=1; subtag='a'; % liver data 2 class (2966 samples): Normal (1277) & HCC (1689)
% tag=1; subtag='b'; % thyroid data 2 class (423 samples): NL (161) & WIFC (262)
% tag=1; subtag='c'; % mesothelioma data 2 class (1080 samples): Benign (590) & Malignant (490)
% tag=1; subtag='d'; % melanoma data 2 class (11542 samples): DN (5189) & MM (6353)

%%
p0=pwd; cd ..
inp_tr=[pwd '/DATA/data' num2str(tag) subtag '/rcdt'];
inp_te=[pwd '/DATA/data' num2str(tag) subtag '/rcdt'];

load([inp_tr '/Rcdt']);
data_tr=u; I0_tr=I0; label_tr=label;

load([inp_te '/Rcdt']);
data_te=u; label_te=label;

cd(p0);

%%
ret_dir_num=4; % 1,2,3....as much directions (features) as you need.

%%
tic
[PCA_directions,projPCA_tr,viz_pca]=PCA_ModesR(data_tr,I0,ret_dir_num);
projPCA_te=proj_n_plot_PCAr(data_te,PCA_directions,viz_pca);

[PLDA_directions,projPLDA_tr,viz_plda]=PLDA_ModesR(data_tr,label_tr,I0_tr,1e+4,ret_dir_num);
[projPLDA_te,f01,f02]=proj_n_plot_PLDAr(data_te,label_te,PLDA_directions,viz_plda);

% fh3=figure('position',[0 40 650 300]); %Linear interpolation in the Radon-CDT space
% index=14;
% alpha=linspace(1,0,5);
% for i=1:length(alpha)    
%     ualpha=alpha(i)*u(:,:,index);    
%     [Ialpha,Ralpha]=iRCDT(ualpha,I0);    
%     subplot(2,5,i)
%     imshow(Ralpha,[]);axis off;axis image    
%     subplot(2,5,5+i)
%     imshow(Ialpha,[]);axis off;axis image   
% end
% [ax1,h1]=suplabel('RI_\alpha=f_\alpha''RI_1(f_\alpha,\theta),  RI_i:Radon transform of I_i','x');
% [ax2,h2]=suplabel('Linear interpolation in the Radon-CDT space','t');
% set(h1,'fontsize',24)
% set(h2,'fontsize',24)
toc
