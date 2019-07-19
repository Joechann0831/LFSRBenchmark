%%% This code is used to perform PCA+RR+BM to super-resolve the light field
%%% images from HCI datasets
clear;
clc;


addpath('./MATLAB/');
addpath('./MATLAB/light_field/');
addpath('./MATLAB/superresolution/');

datasets = {'scene1','scene2','scene3','scene4','scene5','scene6','scene7','scene8','scene9','scene10',...
    'scene11','scene12',};

scales = [2,3];
for mf = 1:2
    scale = scales(mf);
    disp(['----Scaling factor is ',num2str(scale)]);
    for ds = 1:2 % ds = 1,bicubic; ds = 2, gaussian
        if ds == 1
            disp('-------Downsampling method: bicubic');
            ds_method = 'bicubic';
        else
            ds_method = 'gaussian';
            disp('-------Downsampling method: gaussian');
        end
        for i = 1:12
            %% 1, read the HCI data

            lf_name = datasets{1,i};
            disp(['------------Processing ',lf_name]);

            load(['../data/EPFL/',lf_name,'.mat']);

            HR_LF = permute(lf_data,[3,4,5,1,2]);
            
            if ds == 1
                LR_LF = lf_downsample_bicubic(lf_data, scale);
            else 
                LR_LF = lf_downsample_gauss(lf_data, scale);
            end
            
            LR_LF = permute(LR_LF,[3,4,5,1,2]);

            %% 2, bm_pca_rr for the LR_LF

            SR_LF = bm_pca_rr_parent(LR_LF,scale);

            % split the central view to evaluate the results
            Ic_GT = squeeze(HR_LF(:,:,:,5,5));
            Ic_SR = squeeze(SR_LF(:,:,:,5,5));
            Ic_LR = squeeze(LR_LF(:,:,:,5,5));
            % transform the RGB images to gray scale ones
            Ic_GT = rgb2ycbcr(Ic_GT);
            Ic_GT = squeeze(Ic_GT(:,:,1));
            Ic_GT = modcrop(Ic_GT,scale);

            Ic_SR = rgb2ycbcr(Ic_SR);
            Ic_SR = squeeze(Ic_SR(:,:,1));


            Ic_LR = rgb2ycbcr(Ic_LR);
            Ic_LR = squeeze(Ic_LR(:,:,1));

            % calculate the PSNR
            PSNR_SR = psnr(Ic_GT,Ic_SR);
            PSNR_LR = psnr(Ic_GT,Ic_LR);

            % display the PSNR
            disp(['The results of image ',lf_name,' :']);
            disp(['The PSNR of bicubic interpolation is :',num2str(PSNR_LR)]);
            disp(['The PSNR of PCA+RR+BM is :',num2str(PSNR_SR)]);

            % write the image to .mat
            savedir = ['./Results/EPFL/scale',num2str(scale),'/',ds_method,'/'];
            if exist(savedir,'dir') == 0
                mkdir(savedir);
            end
            save([savedir,lf_name,'.mat'],'SR_LF','LR_LF','PSNR_LR','PSNR_SR');
        end
    end
end