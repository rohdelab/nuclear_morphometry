function [ out ] = mean_dev_probe(data_cell,label_cell,cPLDA_directions)
SD_spread=2;
for a=1:length(data_cell)
    label=label_cell{a};
    temp=data_cell{a};
    cPLDA_projections=cPLDA_directions'*temp;
    
    class=sort(unique(label)); pt=15;
    for which=1:size(cPLDA_directions,2)
        for i=1:length(class)
            std_PLDA=SD_spread*std(cPLDA_projections(which,:));
            lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
            [n] =hist(cPLDA_projections(which,label==class(i)),lambdaPLDA');
            XWC(:,i)=n/sum(n);
            VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        end
        mean_values=sum(VWC);
        dm(which,a)=diff(mean_values);
    end
end
s=sign(dm');
mata=[sum(s==-1);sum(s==1)];
s_mata=sum(mata); mx_mata=max(mata);

out=mx_mata./s_mata; out=out(:);
end

