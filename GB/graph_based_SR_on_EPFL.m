%%% This script is used to perform GB-based light field SR methods on HCI
%%% dataset

clear;
clc;

addpath('../functions/');

datasets = {'scene1','scene2','scene3','scene4','scene5','scene6'};

for ind = 1:6
    
    %% 1, read the HCI data
    lf_name = datasets{1,ind};
    scale = 2;
    
    load(['../data/EPFL/',lf_name,'.mat']);
    
    HR_LF = permute(lf_data, [3,4,5,1,2]);
    
    LR_LF = lf_downsample_gauss(lf_data, scale);
    LR_LF = permute(LR_LF, [3,4,5,1,2]);
    
    %% 2, GB-based SR for the LR_LF
    dispMax = 6;
    sub_patch_size = 70;
    tic;
    SR_LF = graph_based_SR(LR_LF, scale, size(HR_LF),dispMax, sub_patch_size);
    elapsed_time = toc;
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
    Ic_LR = imresize(Ic_LR,scale);

    % calculate the PSNR
    PSNR_SR = psnr(Ic_GT,Ic_SR);
    PSNR_LR = psnr(Ic_GT,Ic_LR);
%     save(['HCI_results_gauss/Results_',lf_name,'.mat'],'PSNR_SR','PSNR_LR','SR_LF');
    
    save(['Results/scale',num2str(scale),'/gaussian/EPFL/',lf_name,'.mat'],'PSNR_SR','PSNR_LR','SR_LF','elapsed_time');

    % display the PSNR
    disp(['The results of image ',lf_name,' :']);
    disp(['The PSNR of bicubic interpolation is :',num2str(PSNR_LR)]);
    disp(['The PSNR of GB is :',num2str(PSNR_SR)]);
    
end