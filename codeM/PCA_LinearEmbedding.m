function [Vecmode,eigval]=PCA_LinearEmbedding(Psi);
%% This code generates PCA modes for the linear embedding with 
% respect to the template image. U(:,,i) and V(:,:,i) are the
% translation functions warping image i to the template image.
[nfeature,Nd]=size(Psi);
if nfeature>Nd
    T=(1/Nd)*(Psi'*Psi);
    [eigvec,eigvalue]=eig(T);
    eigvalue(isnan(eigvalue))=0;
    eigvalue(abs(eigvalue)==inf)=1e10;
    eigval = diag(eigvalue);
    [eigval I] = sort((eigval),'descend');    
    Vec=real(eigvec(:,I));    
    %% Calculate modes
    for i=1:size(Vec,2)
        Vecmode(:,i)=(1/sqrt(Nd*eigval(i)))*Psi*Vec(:,i);
    end
else
    ST=(1/Nd)*(Psi*Psi');
    [eigvec,eigvalue]=eig(ST);
    eigvalue(isnan(eigvalue))=0;
    eigvalue(abs(eigvalue)==inf)=1e10;
    eigval = diag(eigvalue);
    [eigval I] = sort((eigval),'descend');
    I(eigval==0)=0;
    Vecmode=real(eigvec(:,I(I~=0)));    
end    
