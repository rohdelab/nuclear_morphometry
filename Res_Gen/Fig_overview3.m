clc
clear all
close all

%%
x = linspace(-4,2,12); sig=1;
gaussian = @(x) (1/sqrt((2*pi*sig*sig))*exp(-x.^2/(2*sig)));
y=gaussian(x);

CL={[0.4 0.65 0.2],[0.9 0.7 0.12],[0.5 0.5 .5]};

%%
f1=figure('Position',[0 0 400 120]); movegui(f1,'north'); hold on
bar(x,y,'facecolor',CL{2})
plot([-1 -1],[0 .5],':','linewidth',2,'color',CL{3})
axis off
