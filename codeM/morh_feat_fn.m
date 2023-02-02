function [tf] = morh_feat_fn(xx)
xx=mat2gray(xx);

%%
x=double(xx>0); x=imfill(x,'holes');
sx = regionprops(x,'area','MajorAxisLength','MinorAxisLength','Eccentricity',...
    'Orientation','ConvexArea','FilledArea','EulerNumber',...
    'EquivDiameter','Solidity','Extent','Perimeter');
sx13=entropy(xx)*100;
sx14=mean(xx(:))*100;
sx15=std(xx(:))*100;
sx16=skewness(xx(:));
sx17=kurtosis(xx(:));

glcms=graycomatrix(xx); 
% t=sort(glcms(:),'descend'); glcms(1)=t(2); 
% sx18 = graycoprops(glcms)
sx18 = graycoprops(glcms(2:end,2:end));

sx20=100*psnr(xx,ones(size(xx)));
sx21=bweuler(xx>0.5*mean(xx(find(xx))));

%%
if length(sx)==0
    tf=zeros(21,1);
else
    tf(1)=sx.Area; tf(2)=sx.MajorAxisLength; tf(3)=sx.MinorAxisLength;
    tf(4)=sx.Eccentricity; tf(5)=sx.Orientation; tf(6)=sx.ConvexArea;
    tf(7)=sx.FilledArea;  tf(8)=sx.EulerNumber; tf(9)=sx.EquivDiameter;
    tf(10)=sx.Solidity; tf(11)=sx.Extent; tf(12)=sx.Perimeter;
    tf(13)=sx13; tf(14)=sx14; tf(15)=sx15; tf(16)=sx16; tf(17)=sx17;
    tf(18)=sx18.Contrast*10; tf(19)=sx18.Energy*100; tf(20)=sx20;
    tf(21)=sx21;
    tf=tf(:);

    tf(isnan(tf))=0; tf(isinf(tf))=0;
end
end