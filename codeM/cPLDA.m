function [ Vec,eigval] = cPLDA(DataB,LabelsB,Alpha,ssall,nPLDA)

%%
ST=0;
SWNew=0;
ss=[];
for a=1:length(DataB)
    Data=DataB{a};
    Labels=LabelsB{a};
    
    %%
    [nfeatures,nsamples]=size(Data);
    if ~(exist('nPLDA'))
        if nfeatures<nsamples
            nPLDA=nfeatures;
        else
            nPLDA=nsamples;
        end
    end
    
    %%
    
    %% FEATURE SCALING PART
    Data=ssall*Data;
    
    
    [ST_s,SWNew_s,ss_s] = cPLDA_p1(Data,Labels,Alpha,nPLDA);
    ST=ST+ST_s;
    SWNew=SWNew+SWNew_s;
    ss{a}=ss_s;
end
% ST=ST/length(DataB);
% SWNew=SWNew/length(DataB);

%%
[Vec,eigval] = cPLDA_p2(ST,SWNew,nPLDA,ssall);

end


