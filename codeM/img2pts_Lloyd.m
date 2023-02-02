function [res_P2,res_c2,llerr,Pl,V,var_out] = img2pts_Lloyd(img,Nmasses)

% This function implements the Voronoi diagram (k-mean) for particle
% approximation. 
% Input: 
% 1. img: the raw gray-level image(image_size_X*image_size_Y);
% 2. Nmasses: number of initial particles (scalor)
% Output:
% 1. res_P2: the locations of the particles
% 2. res_c2: the weights of the particles 
% 3. llerr: the error plot recording this image
% 4. Pl: the locations of all pixels
% 5. V: the pixel intensity value
% 6. var_out: the variations for each Voronoi cell
% Author: Wei Wang

stopLloyd = 0.5;
vis = 0; %visualization?

img = double(img);
[ny,nx]=size(img);

% get threshold using Otsus method
img_t = img./max(img(:));
% img_t = img./sum(img(:));
level = graythresh(img_t)*0.22;
BW = double(im2bw(img_t,level));

% se = strel('disk',10);
% BW = imclose(BW,se);
BW = double(bwareaopen(BW, 7));

STATS = regionprops(BW, 'BoundingBox');
w = find(img_t < level);
img(w) = 0;
ind = find(img_t>=level);
% generate random samples in the domain where there is intensity values
sx = round(STATS.BoundingBox(1));
sy = round(STATS.BoundingBox(2));
Nx = round(STATS.BoundingBox(3));
Ny = round(STATS.BoundingBox(4));
if ~exist('Nmasses','var')
    Nmasses = fix(length(ind)/25)+20;
end

output_Index = randSampling(ind,min([Nmasses,length(ind)]));
res_P = fromInd2Coord(output_Index,ny);
res_c = res_P(1,:)*0;
res_c(:) = 1;


% get coordinates inside bounding box
BW = img_t*0;
BW(:,:) = 0;
img_x = img./sum(img(:));
BW(sy:sy+Ny-1,sx:sx+Nx-1) = img_x(sy:sy+Ny-1,sx:sx+Nx-1);
[row,col,V] = find(BW);

V = V';
Pl = [col'; row']; %[x;y]

if (vis>0)
    figure(2);set(gcf,'Color',[1 1 1])
    imagesc(imresize(img,2*size(img)));colormap gray; axis off;truesize
    hold on
    plot(2*res_P(1,:),2*res_P(2,:),'.r');
    hold off
    drawnow
end

if length(ind)<Nmasses
    
    res_P2 = fromInd2Coord(ind,ny);
    nlz = sum(img(ind));
    res_c2 = (img(ind)./nlz);
    var_out = zeros(1,length(ind));    
    llerr = 0;
    
else

% for it=1:Nit
cur = 1;
differ = 1;
while differ>stopLloyd
    % get nearest neighbors
    neighbors_map = Pl(1,:)*0;
    for k=1:size(Pl,2) %one pass through all pixels in the image (bounding box)
        Pk = Pl(:,k);
        BP = repmat(Pk,1,size(res_P,2));
        err = sum(( (BP-res_P).^2 ),1);
        w = find( err == min(err) );
        neighbors_map(k) = w(1);
    end

    % compute center of mass of each cluster
    for k=1:size(res_P,2)
        w = find(neighbors_map == k);
        cx = sum(  V(w).*Pl(1,w) )./sum(V(w)+1e-10);
        cy = sum(  V(w).*Pl(2,w) )./sum(V(w)+1e-10);
        tmp_center = [cx;cy];
        dist_cent = L2_distance(Pl(:,w),tmp_center).^2;
        errUB(k) = V(w)*dist_cent;
        res_P(:,k) = tmp_center;         
        res_c(k) = sum( V(w) ); %recompute weithgs
    end
    
    llerr(cur) = sum(errUB);
    
    if (vis>0)
        figure(2);set(gcf,'Color',[1 1 1])
        imagesc(imresize(img,2*size(img)));colormap gray; axis off;truesize
        hold on
        plot(2*res_P(1,:),2*res_P(2,:),'.r');
        hold off
        drawnow
    end

if cur>=4
% differ = (sqrt(llerr(cur))-sqrt(llerr(cur-1)))/(sqrt(llerr(cur-1))-sqrt(llerr(cur-2)));
differ = ((llerr(cur))-(llerr(cur-1)))/((llerr(cur-1))-(llerr(cur-2)));
else
differ = 1;
end
cur = cur+1;
end

% neighbors_map = Pl(1,:)*0;
% for k=1:size(Pl,2) %one pass through all pixels in the image (bounding box)
%     P = Pl(:,k);
%     BP = repmat(P,1,size(res_P,2));
%     err = sum(( (BP-res_P).^2 ),1);
%     w = find( err == min(err) );
%     neighbors_map(k) = w(1);
% end
for k=1:size(res_P,2)
    w = find(neighbors_map == k);
    temp = Pl(:,w)-repmat(res_P(:,k),[1,length(w)]);
    vari(k) = std( diag(temp'*temp) ); % Compute the variance
end
%% remove points with coefficient zero
eps = 1e-10;
w = find(res_c < eps);
res_c(w) = [];
res_P(:,w) = [];
vari(w) = [];
var_out = vari;
var_out(var_out==0)=1e-8;

res_c = res_c./sum(res_c);  

res_P2 = res_P;
res_c2 = vec(res_c);

end


end
