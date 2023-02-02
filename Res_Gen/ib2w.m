function [jm] = ib2w(im,df)
i=mat2gray(im);

% %
if df==1
    i(find(i>0.000000001))=1;
elseif df==2
    i(find(i>0.000000001))=1;
else
    i(find(i>0.1))=1;
end

% %

i=im2bw(imfill(i,'holes'));

if df==1
    se = strel('line',5,5);
    i=imdilate(i,se);
elseif df==2
%     se = strel('line',5,5);
%     i=imdilate(i,se);
else
    se = strel('line',20,20);
    i=imerode(i,se);
end


jm=mat2gray(im); jm(find(i==0))=1; jm=imadjust(jm);
end