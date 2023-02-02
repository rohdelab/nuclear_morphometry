function [PLDA_directions,PLDA_projections,viz_plda] = PLDA_ModesR(u,label,I0,alpha,ret_dir_num)
SD_spread=2;
[M,N,K]=size(u);
Psi=zeros(M*N,K);
for k=1:K
    Psi(:,k)=vec(u(:,:,k));
end
Psimean=mean(Psi')';
Psi=Psi-repmat(Psimean,1,K);

alt=2; % 1 or 2
switch alt
    case 1
%         alternate 1 of 2 (upto next newline space)% 1 was original, 2 is borrowed
        [U,S,V]=svd(Psi,'econ');
        PsiPCA=U'*Psi;
        [Vec,~] = PLDA(PsiPCA,label,alpha);
        Vec=real(Vec);
        PLDA_dir=U*Vec;
    case 2
%         alternate 2 of 2 (upto next newline space)% 1 was original, 2 is borrowed
        [EVc,~]=PCA_LinearEmbedding(Psi); EVc=real(EVc);
%         alpha=[];
%         if isempty(alpha)
%             Curveoption.low=.01;
%             Curveoption.high=10;
%             Curveoption.step=0.01;
%             Curveoption.nPLDA=3;
%             alpha=Calculate_Alpha(Psi,label,EVc,Curveoption);
%         end
%         alpha
        alpha=4;%size(Psi,2);
        [Vec,~]=PLDA_in_PCAspace(Psi,label,alpha,size(Psi,2),EVc);
        Vec=real(Vec);
        PLDA_dir=Vec;
    otherwise
        disp('ERROR!!');
end

%%
dir_num=ret_dir_num;% min(size(PLDA_dir,2),2);
PLDA_directions=PLDA_dir(:,1:ret_dir_num);

PLDA_proj=Psi'*PLDA_directions;

for i=1:dir_num
    lambda=linspace(-SD_spread*std(PLDA_proj(:,i)),SD_spread*std(PLDA_proj(:,i)),7);
    big=[];
    for j=1:7
        utemp=reshape(Psimean+lambda(j)*PLDA_directions(:,i),M,N);
        big=[big,mat2gray(iRCDT(utemp,I0))];
    end
    viz_plda(:,:,i)=big;
end
PLDA_projections=PLDA_proj';

%%
class=sort(unique(label));  pt=15;
for which=1:size(PLDA_directions,2)
    for i=1:length(class)
        std_PLDA=SD_spread*std(PLDA_projections(which,:));
        lambdaPLDA=linspace(-std_PLDA,std_PLDA,pt);
        [n] =hist(PLDA_projections(which,label==class(i)),lambdaPLDA');
        XWC(:,i)=n/sum(n);
        VWC(:,i)=XWC(:,i).*[1:length(XWC(:,i))]';
    end
    mean_values=sum(VWC);
    dm(which)=diff(mean_values);
    if dm(which)>0
    else
        temp=PLDA_directions(:,which);
        PLDA_directions(:,which)=-temp;
        temp=viz_plda(:,:,which);
        viz_plda(:,:,which)=fliplr(temp);
    end
end
PLDA_projections=PLDA_directions'*Psi;


