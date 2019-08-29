%%% This script is used to test LFCNN on the whole light field with angular
%%% resolution of 9*9

clear;
clc;

addpath('/home/joechan/caffe/matlab');
addpath('../../functions/');

% Iters = 3000:3000:165000;
% Iters = 129000;

lf_dim = 9;
scale = 2;
ds_method = 'bicubic';
% ds_method = 'gaussian';
modes = 'EPFL1';
dataset = 'scene1';
for iter = 1:length(Iters)
    iter_num = Iters(iter);
    disp(num2str(iter_num));

    %% set devices and load the network

    caffe.set_mode_gpu();
    gpu_id = 0;
    caffe.set_device(gpu_id);

    net_config = '../training_configs/deploy_LFVDSR_net.prototxt';
    net_weights = ['../Models/EPFL1_scale2_bicubic.caffemodel'];
    net = caffe.Net(net_config,net_weights,'test');

    %% load the light fields and test the network, save the result HR light fields

%     disp(['scale is ',num2str(scale),', while downsampling method is ',ds_method]);
%     disp(['Dataset is ',dataset,' with mode ',modes]);

    load(['../../data/EPFL/',dataset,'.mat']);

    % downsample the input light fields and take out Y channel
    seed_view = squeeze(lf_data(1,1,:,:,:));
    seed_view = modcrop(seed_view,scale);
    [H,W,~] = size(seed_view);

    lf_input = zeros([lf_dim,lf_dim,H,W],'uint8');
    lf_gray = lf_rgb2gray(lf_data);
    if strcmpi(ds_method,'bicubic')
        lf_color_ds = lf_downsample_bicubic(lf_data,scale);
    elseif strcmpi(ds_method,'gaussian')
        lf_color_ds = lf_downsample_gauss(lf_data,scale);
    else
        disp('Wrong downsampling kernel name!');
    end
    lf_gray_ds = lf_rgb2gray(lf_color_ds);


    for u = 1:lf_dim
        for v = 1:lf_dim
            % get the downsampled image
            sub_img = squeeze(lf_gray_ds(u,v,:,:));

            % get the input images
            lf_input(u,v,:,:) = imresize(sub_img,scale,'bicubic');
        end
    end

    lf_input = im2double(lf_input);
    LF_SR = zeros(size(lf_input));

    coords = [1,2,2,3,4,5,6,7,8,9];


    for i = 1:1:5
        for j = 1:5
            cnn_input = zeros(H,W,4,1);
            input1 = squeeze(lf_input(coords((i-1)*2+1),coords((j-1)*2+1),:,:));
            input2 = squeeze(lf_input(coords((i-1)*2+1),coords((j-1)*2+2),:,:));
            input3 = squeeze(lf_input(coords((i-1)*2+2),coords((j-1)*2+1),:,:));
            input4 = squeeze(lf_input(coords((i-1)*2+2),coords((j-1)*2+2),:,:));

            cnn_input(:,:,1,:) = input1;
            cnn_input(:,:,2,:) = input2;
            cnn_input(:,:,3,:) = input3;
            cnn_input(:,:,4,:) = input4;

            cnn_input = permute(cnn_input,[2 1 3 4]);
            net.blobs('data').reshape([W H 4 1]);
            net.blobs('data').set_data(cnn_input);
            net.forward_prefilled();

            lf_output = net.blobs('res').get_data();
            lf_output = double(lf_output);
            lf_output = permute(lf_output,[4 3 2 1]);

            LF_SR(coords((i-1)*2+1),coords((j-1)*2+1),:,:) = lf_output(1,1,:,:);
            LF_SR(coords((i-1)*2+1),coords((j-1)*2+2),:,:) = lf_output(1,2,:,:);
            LF_SR(coords((i-1)*2+2),coords((j-1)*2+1),:,:) = lf_output(1,3,:,:);
            LF_SR(coords((i-1)*2+2),coords((j-1)*2+2),:,:) = lf_output(1,4,:,:);
        end
    end
    LF_SR = im2uint8(LF_SR);
    lf_input = im2uint8(lf_input);

    [~,~,~,~,bic_PSNR,~] = PSNR_SSIM(lf_gray,lf_input,scale);
    [ave_PSNR,var_PSNR,ave_SSIM,var_SSIM,cv_PSNR,cv_SSIM] = PSNR_SSIM(lf_gray,LF_SR,scale);
    disp(['PSNR of central view is ',num2str(cv_PSNR),' while bicubic is ',num2str(bic_PSNR)]);
end

savedir = ['./Results/scale',num2str(scale),'/',ds_method,'/EPFL/'];
if exist(savedir,'dir') == 0
    mkdir(savedir);
end
save([savedir,dataset,'.mat'],'LF_SR','cv_PSNR','cv_SSIM');
