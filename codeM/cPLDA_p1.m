function [ST,SWNew,ss] = cPLDA_p1(Data,Labels,Alpha,nPLDA)

[nfeatures,nsamples]=size(Data);
if ~(exist('nPLDA'))
    if nfeatures<nsamples
        nPLDA=nfeatures;
    else
        nPLDA=nsamples;
    end
end

%%
x=mean(Data')';
nclasses=unique(Labels);
SB=zeros(nfeatures,nfeatures);
SW=zeros(nfeatures,nfeatures);
for i=1:length(nclasses)
    Class{i}=Data(:,Labels==nclasses(i));
    mu(:,i)=mean((Class{i})')';
    Sb{i}=(mu(:,i)-x)*(mu(:,i)-x)';
    Sw{i}=zeros(nfeatures,nfeatures);
    ClassiM=repmat(mu(:,i),1,size(Class{i},2));
    Sw{i}=(Class{i}-ClassiM)*(Class{i}-ClassiM)';    
    SB=SB+(1/nsamples)*size(Class{i},2)*Sb{i};
    SW=SW+(1/nsamples)*Sw{i};
end
ST=SB+SW;
SWNew=SW+(Alpha*nsamples)*eye(size(SW));
ss=0;

end


