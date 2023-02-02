function [Pl,P] = particleApproximation(img_array_out, Nmasses,paral)
% This function implements the particle approximation procedure, as
% described in the paper "Analyzing images with linear optimal
% transportation". The basic idea is to first use a K-mean method to
% initially approximate images by particles. Then we either increase the
% number of particles when the error is big (by splitting), or reduce the
% number of particles when the error is small (by particle merging).

% Input:
% 1. img_array_out: the raw gray-level image matrix (3d) with each
% dimension (image_size_X, image_size_Y, N) (N: number_of_image);
% 2. Nmasses: number of initial particles (scalor)
% Output:
% 1. Pl: the locations of the particles (2*N)
% 2. P: the weights of the particles (1*N)
% 3. llerr: the error plot recording in steps for each images (N cells)
% 4. absErr: the absolute error for each approximation (1*N)
% Author: Wei Wang

%% For each image, do the particle approximation by Voronoi diagram (k-mean)
[Ny,Nx,numImg] = size(img_array_out);

if paral
    parfor i = 1:numImg
        display(['Approximating particles for image #',num2str(i)])
        [Pl{i},P{i},llerr{i},Pl_org{i},P_org{i},var_out{i}] = img2pts_Lloyd(img_array_out(:,:,i),Nmasses);
    end
else
    for i = 1:numImg
        display(['Approximating particles for image #',num2str(i)])
        [Pl{i},P{i},llerr{i},Pl_org{i},P_org{i},var_out{i}] = img2pts_Lloyd(img_array_out(:,:,i),Nmasses);
    end
end

