%%% This script is used to test VDSR on the whole light field with angular
%%% resolution of 9*9

clear;
clc;

addpath('/home/joechan/caffe/matlab');
% datasets = {'buddha','buddha2','coneHead','medieval','monasRoom','papillon','rx_elephant','rx_watch','statue','stillLife'};
datasets = {'scene1','scene2','scene3','scene4','scene5','scene6','scene7','scene8','scene9','scene10','scene11','scene12'};
lf_dim = 9;
scale = 3;
ds_method = 'gaussian';

%% set devices and load the network

caffe.set_mode_gpu();
gpu_id = 0;
caffe.set_device(gpu_id);

net_config = '../training_config/VDSR_net_deploy.prototxt';
net_weights = '../Models/_iter_155000.caffemodel';
net = caffe.Net(net_config,net_weights,'test');

%% load the light fields and test the network, save the result HR light fields
for ind = 1:12
    dataset = datasets{ind};
    disp(['Processing light field ',dataset]);
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
            sub_img = squeeze(lf_gray_ds(u,v,:,:));
            lf_input(u,v,:,:) = imresize(sub_img,scale,'bicubic');
        end
    end

    lf_input = im2double(lf_input);
    LF_SR = zeros(size(lf_input));

    for i = 1:1:lf_dim
        for j = 1:lf_dim
            cnn_input = zeros(H,W,1,1);
            cnn_input(:,:,1,1) = squeeze(lf_input(i,j,:,:));

            cnn_input = permute(cnn_input,[2 1 3 4]);
            net.blobs('data').reshape([W H 1 1]);
            net.blobs('data').set_data(cnn_input);
            net.forward_prefilled();

            sub_SR = net.blobs('sum').get_data();
            sub_SR = double(sub_SR);
            sub_SR = permute(sub_SR,[4 3 2 1]);

            LF_SR(i,j,:,:) = sub_SR;
        end
    end
    LF_SR = im2uint8(LF_SR);
    lf_input = im2uint8(lf_input);

    [~,~,~,~,bic_PSNR,~] = PSNR_SSIM(lf_gray,lf_input,scale);
    [ave_PSNR,var_PSNR,ave_SSIM,var_SSIM,cv_PSNR,cv_SSIM] = PSNR_SSIM(lf_gray,LF_SR,scale);
    disp(['PSNR of central view is ',num2str(cv_PSNR),' while bicubic is ',num2str(bic_PSNR)]);
    
    savedir = ['./VDSR_Results/scale',num2str(scale),'/',ds_method,'/EPFL/'];
    if exist(savedir,'dir')==0
        mkdir(savedir);
    end
    save([savedir,dataset,'.mat'],'LF_SR','cv_PSNR','cv_SSIM','lf_input');
end