%%% This script is used to implement the projection based SR method with
%%% the estimated disparity provided by Tingchun Wang's occlusion-aware
%%% disparity estimation 
%%% together with the intensity consistency checking instead of depth consistency
%%% and the back-projection refinement.

clear;
clc;
addpath('../functions/');

datasets = {'scene1','scene2','scene3','scene4','scene5','scene6','scene7','scene8','scene9','scene10','scene11','scene12'};
thresh = 0.10; % threshold parameter for intensity consistency checking
PSNRs = zeros([2,12]);
SSIMs = zeros([2,12]);
for ind_thre = 1:1
    disp(['Threshold is ',num2str(thresh)]);
    for ds_method = 1:2
        if ds_method == 1
            ds_meth = 'bicubic';
            ds_flag = 1;
        else
            ds_meth = 'gaussian';
            ds_flag = 0;
        end
        disp(['------ Scale is 2 while downsampling method is ',ds_meth]);
        for data_ind = 1:12
        %% Input the light field data with only grayscale

            lf_name = datasets{data_ind};
            disp(['---------- Input light field ',lf_name]);
            data_file_name = ['../data/EPFL/',lf_name,'.mat'];
            load(data_file_name);

            scale = 2;
            lf_gray = lf_rgb2gray(lf_data);
            if ds_method == 1
                lf_color_ds = lf_downsample_bicubic(lf_data,scale);
            else
                lf_color_ds = lf_downsample_gauss(lf_data,scale);
            end
            lf_gray_ds = lf_rgb2gray(lf_color_ds);

            [U,V,X,Y] = size(lf_gray_ds);

            %% Estimate the depth of the central view

            disp('Depth estimation start--------------------------');
            disp(' ');
            tic
            lf_data_for_depth = permute(lf_color_ds,[3 4 5 1 2]);
            data_input = im2double(lf_data_for_depth(:, :, :, :, end:-1:1));
            depth_max = 2;
            disparity = computeDepth(data_input , depth_max);
            %% Super resolution by projection method

            disp('Projection-based projection start-----------------');
            sr_pro = sr_projection_scale2(lf_gray_ds , disparity, ds_flag, thresh);
            complete_time = toc;
            ori_ds = squeeze(lf_gray_ds(5 , 5 , : , :));% the downsampled central view
            inter_img = imresize(ori_ds , scale);% bicubic upsample
            central_view = squeeze(lf_gray(5 , 5 , 1:X*scale , 1:Y*scale));% ground truth central view

            psnr_pro = psnr(central_view , sr_pro);
            ssim_pro = ssim(central_view , sr_pro);
            psnr_bic = psnr(central_view , inter_img);
            disp('Super resolution with projection method is done.');
            disp(strcat('PSNR of projection is :',num2str(psnr_pro),',while PSNR of bicubic is :',num2str(psnr_bic)));
            disp(' ');
            
            savedir = ['./Results/EPFL/',ds_meth,'/'];
            if exist(savedir,'dir') == 0
                mkdir(savedir);
            end
            save([savedir,lf_name,'.mat'],'sr_pro','inter_img','central_view');
            PSNRs(ds_method,data_ind) = psnr_pro;
            SSIMs(ds_method,data_ind) = ssim_pro;
        end
    end
end