function img_array_centered = init_center(img_array)
% Center and crop to eliminate the translation
[ny,nx,nz] = size(img_array);
mx = round(nx/2);
my = round(ny/2);

for k=1:nz
    img = single(img_array(:,:,k));
    [cx2,cy2] = center_of_mass2d(img);
    trans.theta = 0;
    trans.tx = 0;
    trans.ty = 0;
    trans.scx = 1;
    trans.scy = 1; 
    trans.sy = 0.0;
    trans.tx = round(cx2 - mx);
    trans.ty = round(cy2 - my);
    trans.cx = 0;
    trans.cy = 0;
    result_b = apply_trans2d(img,trans,1);
    res = result_b(my-round(ny/2-1):my+floor(ny/2),mx-round(nx/2-1):mx+floor(ny/2));
    img_array_centered(:,:,k) = res;
    
%    figure(1);imagesc(img);colormap gray;truesize
%    figure(2);imagesc(result_b);colormap gray;truesize
%    figure(3);imagesc(res);colormap gray;truesize
    %pause
end






function [cx,cy] = center_of_mass2d(img_input)
img_t = img_input./max(img_input(:));
level = graythresh(img_t)*0.8;
BW = double(im2bw(img_t,level));

img = img_input.*BW;
[M,N] = size(img);
[X,Y] = meshgrid(1:N,1:M);

cx = sum(sum(X.*img))/sum(sum(img));
cy = sum(sum(Y.*img))/sum(sum(img));
end



function result = apply_trans2d(img,trans,degree)

%
% function result = apply_trans2d(img,trans,degree)
% NOTE: if Bspline of degrees 2 or higher are used, image must be of size
% 2^N*2^N.
% Author: G.K. Rohde, 10/09/06

[M,N] = size(img);
[X,Y] = build_trans2d(M,N,trans);

if nargin < 3, 
	degree = 1; 
end

% Interpolate
if (degree < 2)	
	if (degree == 0)
		result = interp2(img,X,Y,'nearest',0);
	else 
		result = interp2(img,X,Y,'linear',0);
	end	
else	
	xcoord = reshape(X,M*N,1);
	ycoord = reshape(Y,M*N,1);
	R = my_interpolate_bspline(img,xcoord,ycoord,degree);
	result = reshape(R,M,N);	
end
end



function [X,Y] = build_trans2d(M,N,transf)

%
% function [X,Y] = build_trans2d(M,N,transformation)
%

[Xt,Yt] = meshgrid(1:N,1:M);
a = transf.theta; %should be in radians
R = [cos(a),sin(a);-sin(a),cos(a)];
S = [transf.scx,transf.sy;0,transf.scy];
M2 = S*R;

%cx = N/2;
%cy = M/2;
cx = transf.cx;
cy = transf.cy;

Xt = Xt-cx;
Yt = Yt-cy;

X = M2(1,1)*Xt + M2(1,2)*Yt + transf.tx + cx;
Y = M2(2,1)*Xt + M2(2,2)*Yt + transf.ty + cy;
end






end