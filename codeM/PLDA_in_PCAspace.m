function [P_Vec,sigma]=PLDA_in_PCAspace(Psi,labels,Alpha,nPLDA,VecPCA)
% This function calculates the PLDA in PCA space
% flag determines whether to save the results or not        


[nfeature,Nd]=size(Psi);
Nvec=size(VecPCA,2);
Psi_PCA=VecPCA'*Psi;
class=unique(labels);
%%
[PLDA_directions,var]=PLDA(Psi_PCA,labels,Alpha,nPLDA);
PLDA2=PLDA_directions'*(Psi_PCA);
sigma=sqrt(var)*Nd;
P_Vec=VecPCA*PLDA_directions;
