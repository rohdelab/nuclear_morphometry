clc
clear all
close all

%% TRANSFORMS & PROCESSINGS
% r6_lotp_END % need to clear up some parts, saved lotp are okay for now. match to them
              % in ot_scaling_linear.m lines 70-79 

%% PARAMETER SELECT
TagVec = 1; % 4; %-2; %1:4; % Choose 1 data, 2 data, or all data: TagVec=0; [1,3]; 0:4.
Ndir=10; Nfold=2; Ndir01=5;

reg_str='TOF'; % 'TOF' 'TSOF' 'TSOF_post'
reg_str_both=1; %0; % use both 'TOF' and 'TSOF' if set = 1

which_axisA=[1 2]; which_axisC=[1:4];  which_axisD=1;
pca_or_plda_in_01=1; % 1 = PCA, 0 = PLDA
SD_spread=3;

save(['../DATA/METADATA/params'],'Ndir','Ndir01','Nfold','reg_str','which_axisA','which_axisC','which_axisD','TagVec','reg_str_both','pca_or_plda_in_01');
save(['../DATA/METADATA/params_inside'],'SD_spread');


%% TRAIN
main_train % A, B, C
main_trainD
main_train01 % A, B, C % run twice with pca_or_plda_in_01 = 0 and 1 

%% TEST
main_testA
main_testB
main_testC
main_testD

main_test01A % run twice with pca_or_plda_in_01 = 0 and 1 
main_test01B % run twice with pca_or_plda_in_01 = 0 and 1 
main_test01C % run twice with pca_or_plda_in_01 = 0 and 1 

