function [ out ] = ttest2_probe(data_cell,label_cell,cPLDA_directions)
h_big=[]; p_big=[]; H_big=[]; P_big=[];
for a=1:length(data_cell)
    label=label_cell{a};
    temp=data_cell{a};
    cPLDA_projections=cPLDA_directions'*temp;
    
    class=sort(unique(label)); x1=[]; x2=[];
    for which=1:size(cPLDA_directions,2)
        x1(:,which)=cPLDA_projections(which,label==class(1));
        x2(:,which)=cPLDA_projections(which,label==class(2));
        
        [h_tmp,~] = ttest2(x1(:,which),x2(:,which));
        
        h_big(a,which)=h_tmp;
    end
end
s=h_big;
mata=[sum(s==0);sum(s==1)];
s_mata=sum(mata);

out=mata(2,:)./s_mata; out=out(:);
end

