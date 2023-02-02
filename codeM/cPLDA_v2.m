function [ Vec,eigval] = cPLDA_v2(DataB,LabelsB,Alpha,nPLDA)

%%
Vec=zeros(size(DataB{1},1));; eigval=0;
for a=1:length(DataB)
    Data=DataB{a};
    Labels=LabelsB{a};
    clss=unique(Labels);

    %%
    [nfeatures,nsamples]=size(Data);
    if ~(exist('nPLDA'))
        if nfeatures<nsamples
            nPLDA=nfeatures;
        else
            nPLDA=nsamples;
        end
    end

    [v_tmp,e_tmp] = PLDA(Data,Labels,Alpha);

    for b=1:size(v_tmp,2)
        v=v_tmp(:,b);
        proj=v'*Data;
        m1=mean(proj(find(Labels==clss(1)))); m2=mean(proj(find(Labels==clss(2))));
        if m1>m2 % abs(m1)<abs(m2)
            v_tmp(:,b)=-v;
        end
    end

    Vec(:,1:size(v_tmp,2))=Vec(:,1:size(v_tmp,2))+v_tmp;

    for i=1:size(Vec,2)
        Vec(:,i)=Vec(:,i)/sqrt(Vec(:,i)'*Vec(:,i));
    end


    eigval=eigval+e_tmp;
end
eigval=eigval/length(DataB);

%%
% for i=1:size(Vec,2)
%     Vec(:,i)=Vec(:,i)/sqrt(Vec(:,i)'*Vec(:,i));
% end

end


