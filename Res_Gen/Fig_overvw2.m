clc
clear all
close all

%%
x = linspace(-1,3,500); 

%%
alpha=[4,-6,2];
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewedgaussian = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);

d1=skewedgaussian(x-0.7,alpha(1)); distribution1=d1(:)*d1(:)';
d2=skewedgaussian(x-1.0,alpha(2)); distribution2=d2(:)*d2(:)';
d3=skewedgaussian(x-1.0,alpha(3)); distribution3=d3(:)*d3(:)';

%%
% distribution1=flipud(distribution1);
% distribution2=fliplr(flipud((distribution2)));
distribution3=fliplr(distribution3);

distribution1=imrotate(distribution1,45,'bilinear','crop');


C2=[0 0.4470 0.7410]; C1=[0.8500 0.3250 0.0980]; 
% C2=[0.65 0.65 .75]; C1=[.85 0.7 0.7]; 
C3=[0.1 0.5 0.1];

Ncl=3; scN=50; msz=300;
falp=.3;

%%
f1=figure('position',[0 0 900 800]); movegui(f1,'north')
[c1,m1]=contour(distribution1,Ncl,'--','color',C1,'linewidth',2); hold on; grid on
[c2,m2]=contour(distribution2,Ncl,'-.','color',C2,'linewidth',2); hold on; grid on
[c3,m3]=contour(distribution3,Ncl,':','color',C3,'linewidth',2); 


[x,y]=find(distribution1>.1); 
ind=randsample(1:length(x),scN); x=x(ind); y=y(ind);
scatter(y,x,msz,'markerfacecolor',C1,'markeredgecolor','k','markerfacealpha',falp);

[x,y]=find(distribution2>.1); 
ind=randsample(1:length(x),scN); x=x(ind); y=y(ind);
scatter(y,x,msz,'markerfacecolor',C2,'markeredgecolor','k','markerfacealpha',falp);

[x,y]=find(distribution3>.1); 
ind=randsample(1:length(x),scN); x=x(ind); y=y(ind);
scatter(y,x,msz,'markerfacecolor',C3,'markeredgecolor','k','markerfacealpha',falp);


axis square; axis off
