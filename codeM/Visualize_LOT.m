function [I1] = Visualize_LOT(Data,Intensity,Nx,Ny,scale,crop)
%This function take a mass point function and generate an image 
%from it haphazardly.
%
%Input:
%      Data= this is a vector of the coordinates
%      [x(1),...,x(n),y(1),...,y(n)]^T.
%      Intensity= this is the vector of the intensity of each coordinate
%      [I(1),...,I(n)]^T.
%      [Nx, Ny]=  Size of the output image
%Output: 
%      I1= Generated Image
%
%Author: Soheil Kolouri
%Date:   09/26/2012
if ~(exist('scale'))
    scale=1;
end 
if ~(exist('crop'))
    crop=1;
end
NG=35; % The size of the Gaussian window 
I1=zeros(scale*Nx,scale*Ny); 
loc=scale*round(reshape(Data,length(Data)/2,2));
linearind=sub2ind(size(I1),loc(:,2),loc(:,1));
I1(linearind)=Intensity; 
h = gausswin(NG*scale)*gausswin(NG*scale)';
h = h./sum(h(:));
I1 =    mat2gray(imfilter(I1,h,'same'));
end
