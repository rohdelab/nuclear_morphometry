function [MM_all,MM_label_all,HH_all] = patient_hist_mean_feature(projPLDA_inp,lab_inp,pat_lab_inp,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
SD_spread=2; pt=16;
for a=1:N
    if N==1
        proj=projPLDA_inp; lab=lab_inp; pat_lab=pat_lab_inp;
        pat_class=unique(pat_lab_inp);
    else
        proj=projPLDA_inp{a}; lab=lab_inp{a}; pat_lab=pat_lab_inp{a};
        pat_class=unique(pat_lab_inp{a});
    end
    
    clear XWC VWC MM
    HH=[];
    for which=1:size(proj,1)
        clear MM_label
        for i=1:length(pat_class)
            std_PLDA=SD_spread*std(proj(which,:));
            lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
            [n] =hist(proj(which,pat_lab==pat_class(i)),lambdaPLDA');
            XWC(:,which,i)=n/sum(n);
            VWC(:,which,i)=XWC(:,which,i).*[1:length(XWC(:,which,i))]';
            MM(which,i)=sum(VWC(:,which,i));
            MM_label(i)=unique(lab(pat_lab==pat_class(i)));
        end
        HH=[HH;squeeze(XWC(:,which,:))];
    end
    if N==1
        XWC_all=XWC; VWC_all=VWC; MM_all=MM; MM_label_all=MM_label;
        HH_all=HH;
    else
        XWC_all{a}=XWC; VWC_all{a}=VWC; MM_all{a}=MM; MM_label_all{a}=MM_label;
        HH_all{a}=HH;
    end
end

end

