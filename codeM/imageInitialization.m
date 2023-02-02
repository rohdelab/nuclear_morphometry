function img_array_out = imageInitialization(img_array, option)

% Eliminate -- Translation, Scaling, Orientation, Flipping (T,S,O,F)
% option = 7 --> 0111 --> Scaling, Orientation, Flipping (no Translation)
% option = 10 --> 1010 --> Translation, Orientation (no Scaling, no Flipping)
% Modified by M. Shifat-E-Rabbi


if option<16


    if ~exist('option','var') % Default parameters of LDDMM
        flag_tran=1;
        flag_scale=1;
        flag_orient=1;
        flag_group=1;
    else    % User inputed parameters
        flag_tran=bitget(option, 4);
        flag_scale=bitget(option,3);
        flag_orient=bitget(option, 2);
        flag_group=bitget(option, 1);
    end

    % Centering the images such that the output objects are located at image center
    if flag_tran==1
        disp('Centering to eliminate the translation factor. Please wait ...');
        img_array_centered = init_center(img_array);
        %save('GL_centered.mat','img_array_centered')
        img_array = img_array_centered;
    end

    % Rescale the images such that the output objects have the same number of
    % pixels
    if flag_scale==1
        % To eliminate the scaling factor
        %         disp('Rescaling images to make the foreground have same area. Please wait ...');
        %         img_array_scaled =  init_scaled(img_array);

        % % %

        series_use=0;

        if series_use==1
            disp('Rescaling images to make the foreground of images have same area using v2. Please wait ...');
            img_array_scaled =  init_scaled_v2(img_array); %save('GL-scaled.mat','img_array_scaled');
            img_array = img_array_scaled;
        else
            disp('Rescaling images to make the foreground of images have same area using v3. Please wait ...');
            tmp=[];
            parfor sind=1:size(img_array,3)
                tmp(:,:,sind) =  init_scaled_v2(img_array(:,:,sind));
                disp(['scale - ' num2str(sind) ' - ' num2str(size(img_array,3)) '. Perc: ' num2str(100*sind/size(img_array,3)) ' %..']);
            end
            img_array=tmp;
        end



        % % %
    end

    % Rotate the images such that the major axis are aligned verticlely
    if flag_orient==1
        % MA align to eliminate orientation, images are rotated to have the
        % same orientation through a principal axis (Hotteling) transform.
        [ny,nx,nz] = size(img_array);
        disp('Aligning images along one axis to eliminate the orientation. Please wait ...');
        hW = waitbar(0,'Aligning images. Please wait ...');
        for k=1:nz,
            img = img_array(:,:,k);
            img_array_ma(:,:,k) = ma_align2d(img);
            waitbar(k/nz);
        end
        close(hW);
        %save('GL-ma.mat','img_array_ma');
        img_array = img_array_ma;
    end

    if flag_group==1
        % Group registration: images are 'flipped' left to right, and up and
        % down for the group affine registration

        %
        %         disp('Flipping images for group registration. Please wait ...');
        %         img_array_out = group_register2d_rigid(img_array);
        disp('Flipping images fast. Please wait ...');
        img_array_out = fast_flip(img_array);
        %

    else
        img_array_out = img_array;
    end

else
    disp('Option value should vary between 0 and 15');
    img_array_out = [];
end