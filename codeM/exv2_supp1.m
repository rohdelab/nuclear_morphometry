function [l] = exv2_supp1(cPLDA_dir,PsiB_m0,labelB,SD_spread,lam)

Ncp=10; cPLDA_dir=cPLDA_dir(:,1:Ncp);
for a=1:length(PsiB_m0)
    t=PsiB_m0{a};
    prj_tmp_tr{a}=cPLDA_dir'*t;
end

for a=1:length(labelB)
    ltmp=labelB{a}; prj_tmp=prj_tmp_tr{a};
    class=sort(unique(ltmp)); pt=15;
    for which=1:size(cPLDA_dir,2)
        for i=1:length(class)
            std_PLDA=SD_spread*std(prj_tmp(which,:));
            lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
            [n] =hist(prj_tmp(which,ltmp==class(i)),lambdaPLDA');
            XWC(:,i)=n/sum(n);
            VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
        end
        mtmp=sum(VWC);
        dmtmp(which,a)=diff(mtmp);
    end
end

dm2=dmtmp;

%


% dm2(find(abs(sum(sign(dmtmp)'))<length(labelB)),:)=0;
% [m,l]=max(sum(abs(dm2')));

% dm2(find(abs(sum(sign(dmtmp)'))<length(labelB)),:)=0;
% [m,l]=max(abs(dm2(:,find(sum(abs(dm2))==min(sum(abs(dm2)))))));
% if length(l)>1
%     l=l(1);
% end


l=min(find(abs(sum(sign(dmtmp)'))==length(labelB)));

end

