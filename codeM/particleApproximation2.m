function [Pl,P,llerr,absErr] = particleApproximation(img_array_out, Nmasses)
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
for i = 1:numImg
    [Pl{i},P{i},llerr{i},Pl_org{i},P_org{i},var_out{i}] = img2pts_Lloyd(img_array_out(:,:,i),Nmasses);
    [Pl{(i)},P{(i)},error,vari{i}] = particleCheckSplitting(Pl{(i)},P{(i)},Pl_org{(i)},P_org{(i)});
    tmpErr = llerr{(i)};
    tmpErr(end+1:end+length(error)) = error; % Update the error after particle-splitting
    llerr{(i)} = tmpErr;
end

%% For all the data, do particle merging to reduce the number of particles (optional)
% for i = 1:numImg
% merError = 1.05*sqrt(llerr{i}(end));
% [Pl{(i)},P{(i)},error] = particleMerging(Pl{(i)},P{(i)},merError^2,Pl_org{(i)},P_org{(i)});
% tmpErr = llerr{(i)};
% tmpErr(end+1:end+length(error)) = error;
% llerr{(i)} = tmpErr;
% end
% 
% saveName = 'particleAppr.mat';
% save(saveName, 'Pl','P','llerr'); 

%% Estimating the errors for the whole data set
% After estimating the error distribution, we will do particle merging 
% (reduce the number of particles) for images with small errors. And We do
% particle splitting (increase the number of particles) for images with
% big errors.
for i = 1:numImg
    absErr(i) = llerr{i}(end);
end
% sigma is the std dev of the particle approximation error distribution.
% The smaller the value, the stricter we want for the error distribution.
sigma = 1; 
absErr = sqrt(absErr);
mean1 = mean(absErr); % Mean of all the errors
cutf1 = std(absErr); % One std deviation of errors
smallErr = mean1-sigma*cutf1;
bigErr = mean1+sigma*cutf1;

smInd=find(absErr<smallErr);
bgInd=find(absErr>bigErr);


%% For small-error data, do the particle merging.
for i = 1:length(absErr)
    if absErr(i)<smallErr
            [Pl{i},P{i},error] = particleMerging(Pl{i},P{i},smallErr.^2,Pl_org{i},P_org{i});
    elseif absErr(i)>bigErr
        [Pl{i},P{i},error] = particleSplitting(Pl{i},P{i},bigErr.^2,Pl_org{i},P_org{i});%,img50(:,:,bgInd(i)));    
    end
    tmpErr = llerr{i};
    tmpErr = [tmpErr,error];
    llerr{i} = tmpErr;
end

for i = 1:numImg
    absErr(i) = llerr{i}(end);
end
absErr = sqrt(absErr);

