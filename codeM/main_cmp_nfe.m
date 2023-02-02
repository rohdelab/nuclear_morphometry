clc
clear all
close all

%% TRANSFORMS & PROCESSINGS
% r6_lotp_END % need to clear up some parts, saved lotp are okay for now. match to them
              % in ot_scaling_linear.m lines 70-79 

%% PARAMETER SELECT
TagVec = 3; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=10; Nfold=2; Ndir01=5;

reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'

which_axisA=[1 1]; which_axisC=[1:1];  which_axisD=1;
% which_axisA=[1 2]; which_axisC=[1:4];  which_axisD=1;
pca_or_plda_in_01=1; % 1 = PCA, 0 = PLDA
SD_spread=3;

save(['../DATA/METADATA/nfe_params'],'Ndir','Ndir01','Nfold','reg_str','which_axisA','which_axisC','which_axisD','TagVec','pca_or_plda_in_01');
save(['../DATA/METADATA/nfe_params_inside'],'SD_spread');


%% TRAIN
% main_train_cmp_nfe % A, B, C
main_train01_cmp_nfe % A, B, C % run twice with pca_or_plda_in_01 = 0 and 1 

%% TEST
% main_testA_cmp_nfe
main_test01A_cmp_nfe % run twice with pca_or_plda_in_01 = 0 and 1 
