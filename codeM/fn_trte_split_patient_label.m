function [label_tr,label_te]=fn_trte_split_patient_label(label,ind,fold_tr,fold_te)
tmp2=[];
for a=1:length(fold_tr)
    tmp=ind{fold_tr(a)}; tmp=tmp(:);
    tmp2=[tmp2;tmp];
end
ind_tr=tmp2;
label_tr=label(ind_tr);

tmp2=[];
for a=1:length(fold_te)
    tmp=ind{fold_te(a)}; tmp=tmp(:);
    tmp2=[tmp2;tmp];
end
ind_te=tmp2;
label_te=label(ind_te);
end

