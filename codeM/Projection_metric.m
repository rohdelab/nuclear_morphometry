function [d] = Projection_metric( A,B )
%This function calculates the projection metric based on the 
%principal angles between the columns of matrix A and B
%Input:
%      A: Orthonormal matrix N*M 
%      B: Orthonormal matrix N*M
%Output: 
%      d: Distance between column spaces of A and B
%Author: 
%      Soheil Kolouri
%Date: 
%      09/27/2012
[N,M]=size(A);
A=orth(A);
B=orth(B);
[~,S,~]= svd(A'*B);
costheta=diag(S);
d=sqrt(M-sum(costheta.^2));
end

