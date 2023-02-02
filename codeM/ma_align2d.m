function [result] = ma_align2d(img_in);
% Aligning each image by rotating the image such that the major axis is
% verticle. (Using PCA method)

%% Threshold the Image to get the Binary Code
level = graythresh(img_in)*0.8;
img = double(im2bw(img_in,level));

%% PCA alignment
[ny,nx] = size(img); 
[X,Y] = find(img); %Find the coordinates of nonzero elements of the image
X = X-mean(X);     %Center the image
Y = Y-mean(Y);     %Center the image     

% [COEFF,SCORE] = princomp([X,Y]);%Calculate the PCA components of X and Y
[COEFF,SCORE] = pca([X,Y]);%Calculate the PCA components of X and Y


%Note that in princomp(X) rows of X correspond to observations and columns to variables.

v = COEFF(:,1); %Get the first eigenvector 

X = [nx/2;nx/2+50*v(2)];
Y = [ny/2;ny/2+50*v(1)];

%% We want the first eigenvector to be aligned with e1=[1;0].
% Therefore we first calculate the angle between the vectors v and e1.
a = acos(v'*[1;0]); %a is the angle

if (v'*[0;1] < 0)   %cos(a) and cos(-a) are equal
    a = -a;
end

%% Generate the 2D Affine registration matrix
trans.theta = 0;
trans.tx = 0;
trans.ty = 0;
trans.scx = 1;
trans.scy = 1;
trans.sy = 0.0;
trans.cx = nx/2;
trans.cy = ny/2;
trans.theta = a;

%% Apply the transformation to the image 
result_b = apply_trans2d(img_in,trans,1);



%% skewness alignment
res = result_b;
[Y,X] = find(result_b);
X = X;
Y = Y;
skx = skewness(X);
sky = skewness(Y);

if (skx < 0)
    res = fliplr(result_b);
end
if (sky > 0)
    res = flipud(result_b); 
end

result = res;
%%

function result = apply_trans2d(img,trans,degree)
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