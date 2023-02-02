function [xx_tr,label_tr,xx_te,label_te]=fn_trte_split2(xx,label,ind,fold_tr,fold_te)
tmp2=[];
for a=1:length(fold_tr)
    tmp=ind{fold_tr(a)}; tmp=tmp(:);
    tmp2=[tmp2;tmp];
end
ind_tr=tmp2;
xx_tr=xx(:,ind_tr); label_tr=label(ind_tr);

tmp2=[];
for a=1:length(fold_te)
    tmp=ind{fold_te(a)}; tmp=tmp(:);
    tmp2=[tmp2;tmp];
end
ind_te=tmp2;
xx_te=xx(:,ind_te); label_te=label(ind_te);
end

