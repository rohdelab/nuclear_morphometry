function [f,df,u,du,RCDT]=RCDT(I0,I1)
% M Shifat-E-Rabbi

epsilon=1e-6; theta=0:179;

R1=radon(I1,theta);
R0=radon(I0,theta);

f=zeros(size(R1));df=f;u=f;du=f;
for a=1:length(theta)
    J1=R1(:,a)+epsilon;
    J0=R0(:,a)+epsilon;
    [f(:,a),df(:,a),u(:,a),du(:,a)]=CDT(J0,J1);
end

RCDT=u.*sqrt(R0+epsilon);
