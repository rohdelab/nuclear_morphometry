function [Iu,Ru]=iRCDT(u,I1)
% M Shifat-E-Rabbi

epsilon=1e-6; theta=0:179;

[X0,Y0]=meshgrid(1:size(u,2),1:size(u,1));
R1=radon(I1,theta);
[~,uy]=gradient(u);
f=Y0-u; df=1-uy;
Ru=df.*interp2(X0,Y0,R1+epsilon,X0,f);

Ru=inpaint_nans(Ru);

Iu=iradon(Ru-epsilon,theta);
Iu(Iu<0)=0;



