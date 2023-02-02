function [Vec,eigval] = cPLDA_p2(ST,SWNew,Ndir,ss)

[eigvec,eigvalue]=eig(ST,SWNew);
eigvalue(isnan(eigvalue))=0;
eigvalue(abs(eigvalue)==inf)=1e10;
eigval = diag(eigvalue);
[eigval I] = sort((eigval),'descend');
Vec=real(eigvec(:,I(1:Ndir)));

%% FEATURE SCALING INVERSION PART
Vec=inv(ss)*Vec;

%%
for i=1:size(Vec,2)
    Vec(:,i)=Vec(:,i)/sqrt(Vec(:,i)'*Vec(:,i));
end
end


