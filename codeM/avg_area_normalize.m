function img_array_scaled = avg_area_normalize(img_array_centered)
% Rescale the images such that the output objects have the same number of
% pixels

img_array_centered = single(img_array_centered);


%
img_area=[];
for a=1:size(img_array_centered,3)
    temp=img_array_centered(:,:,a)>0;
    img_array_centered_mask(:,:,a)=double(imfill(temp,'holes'));
    img_area(a)=sum(sum(img_array_centered_mask(:,:,a)));
end
avg_img_area=mean(img_area);
%


% To eliminate the scaling factor

%
% pa = sum(sum(img_array_centered(:,:,2)));
perc_area=0.07;
target_avg_img_area=size(img_array_centered,1)*size(img_array_centered,2)*perc_area;
mlt_fct=target_avg_img_area/avg_img_area;
pa_first_image=sum(sum(img_array_centered_mask(:,:,1)))*mlt_fct;
%


trans.theta = 0;
trans.tx = 0;
trans.ty = 0;
trans.scx = 1;
trans.scy = 1;
trans.sy = 0.0;
[ny,nx,nz] = size(img_array_centered);
NDS = 300;
ds = linspace(0,0.5,NDS);
dsp = linspace(0,1,NDS);
trans.cx = nx/2;
trans.cy = ny/2;
hW = waitbar(0,'Rescaling images to make the foreground have same area. Please wait ...');
for k=1:nz,
    img = img_array_centered(:,:,k);
    img_array_scaled(:,:,k) = img;
    
    %
    img = img_array_centered_mask(:,:,k);
    img0 = img_array_centered(:,:,k);
    %
    
    if k==1
        ca =sum(sum(img));
        if (ca < pa_first_image)
            minv = abs(ca-pa_first_image);
            minl = 0;
            for m=1:NDS
                trans.scx = 1-ds(m);
                trans.scy = 1-ds(m);
                result_b = apply_trans2d(img,trans,1);
                cv = sum(sum(result_b));
                if ( abs(cv-pa_first_image) < minv)
                    minv = abs(cv-pa_first_image);
                    minl = ds(m);
                end
            end
            trans.scx = 1-minl;
            trans.scy = 1-minl;
            
            %
            %         result_b = apply_trans2d(img,trans,1);
            result_b = apply_trans2d(img0,trans,1);
            %
            
            img_array_scaled(:,:,k) = result_b;
            
        else
            if (ca > pa_first_image)
                minv = abs(ca-pa_first_image);
                minl = 0;
                for m=1:NDS
                    trans.scx = 1+dsp(m);
                    trans.scy = 1+dsp(m);
                    result_b = apply_trans2d(img,trans,1);
                    cv = sum(sum(result_b));
                    if ( abs(cv-pa_first_image) < minv)
                        minv = abs(cv-pa_first_image);
                        minl = dsp(m);
                    end
                end
                trans.scx = 1+minl;
                trans.scy = 1+minl;
                
                %
                %             result_b = apply_trans2d(img,trans,1);
                result_b = apply_trans2d(img0,trans,1);
                %
                
                img_array_scaled(:,:,k) = result_b;
            end
        end
    else
        result_b = apply_trans2d(img0,trans,1);
        img_array_scaled(:,:,k) = result_b;
    end
    
    
    waitbar(k/nz);
end
close(hW)



    function result = apply_trans2d(img,trans,degree)
        
        
        [M,N] = size(img);
        [X,Y] = build_trans2d(M,N,trans);
        
        if nargin < 3,
            degree = 1;
        end
        
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
        
        [Xt,Yt] = meshgrid(1:N,1:M);
        a = transf.theta; %should be in radians
        R = [cos(a),sin(a);-sin(a),cos(a)];
        S = [transf.scx,transf.sy;0,transf.scy];
        M2 = S*R;
        
        cx = transf.cx;
        cy = transf.cy;
        
        Xt = Xt-cx;
        Yt = Yt-cy;
        
        X = M2(1,1)*Xt + M2(1,2)*Yt + transf.tx + cx;
        Y = M2(2,1)*Xt + M2(2,2)*Yt + transf.ty + cy;
    end



end