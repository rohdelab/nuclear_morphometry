function [ out1,out2 ] = spread_probe(data_cell,label_cell,cPLDA_directions)
h_big=[]; p_big=[]; H_big=[]; P_big=[];
for a=1:length(data_cell)
    label=label_cell{a};
    temp=data_cell{a};
    cPLDA_projections=cPLDA_directions'*temp;
    
    class=sort(unique(label)); x1=[]; x2=[];
    for which=1:size(cPLDA_directions,2)
        x1(:,which)=cPLDA_projections(which,label==class(1));
        x2(:,which)=cPLDA_projections(which,label==class(2));
        
        spread(a,which)=std(cPLDA_projections(which,:));

    end
    
    spread(a,:)=100*spread(a,:)/sum(spread(a,:));
end

out1=spread';
out2=cumsum(spread');
end

