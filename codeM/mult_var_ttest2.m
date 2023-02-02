function [H,P] = mult_var_ttest2(label,cPLDA_projections)
class=sort(unique(label));
if length(class)==2
    x1=[]; x2=[];
    for which=1:size(cPLDA_projections,1)
        x1(:,which)=cPLDA_projections(which,label==class(1));
        x2(:,which)=cPLDA_projections(which,label==class(2));
%         [h_tmp,p_tmp] = ttest2(x1(:,which),x2(:,which));
%         h_big(which)=h_tmp;
%         p_big(which)=p_tmp;    
    end
    szx1=size(x1); szx2=size(x2);
    
    X_HOT=[ones(szx1(1),1) x1;2*ones(szx2(1),1) x2];
    [H,P] = T2Hot2iho_mine(X_HOT);
else
    disp('Skipped multi-variate t-test...')
    H=nan; P=nan;
end
end

