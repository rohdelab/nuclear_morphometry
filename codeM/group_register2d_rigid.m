function img_array_out = group_register2d_rigid(img_array_in)

% Group registration: All the images are flipped up-down and left-right. To
% make sure the L2 norm between the current image and the current template
% are minimized

Num_it = 10;

[ny,nx,nz] = size(img_array_in);
A = cell(nz,1);
trans.theta = 0;
trans.tx = 0;
trans.ty = 0;
trans.scx = 1;
trans.scy = 1;
trans.sy = 0.0;
trans.cx = nx/2;
trans.cy = ny/2;

for k=1:nz
    A{k} = trans;
end

img_array_out = img_array_in;
img_array_temp = img_array_in;

wH = waitbar(0,'Taking group registration. Please wait...');
for k=1:Num_it    
    for i=1:nz

        %update img_array_temp       
        mI = mean(img_array_temp,3);
        
        img = img_array_temp(:,:,i);
        imgc = img;
        imgfr = fliplr(img);
        imgfd = flipud(img);
        imgfrd = fliplr(imgfd);
        
        
        ac = sum(sum( (img-mI).^2 ));
        
        acr = sum(sum( (imgfr-mI).^2 ));
        if (acr < ac)
            imgc = imgfr;
        end
        
        acd = sum(sum( (imgfd-mI).^2 ));
        if (acd < ac)
            imgc = imgfd;
        end
        
        acdr = sum(sum( (imgfrd-mI).^2 ));
        if (acdr < ac)
            imgc = imgfrd;
        end

        img_array_temp(:,:,i) = imgc;
                
    end        
    waitbar(k/Num_it)
end
close(wH)

img_array_out = img_array_temp;


%-------------------------------------------------------------------------
function ar_out = update_img_array_temp(ar_in,A)

[ny,nx,nz] = size(ar_in);

for k=1:nz    
    img = ar_in(:,:,k);    
    ar_out(:,:,k)=apply_trans2d(img,A{k},1);    
end

%-------------------------------------------------------------------------
function [tx,ty,angle] = search_grr(img,trans,mI)

trans_temp = trans;

% search translation x

% search translation y

% searc angle

% Flip


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