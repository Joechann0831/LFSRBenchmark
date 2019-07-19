%%% This script is used to generate the data for LFVDSR
%%% The original code is forked from Hanzhi
%%% We need to change the downsample method and scale. So there are amout
%%% of data needed to generate here, we generate the mat file first and
%%% get the h5 files later
%%% 
clear;
clc;

set_names = {'buddha','buddha2','coneHead','medieval','monasRoom','papillon','rx_elephant','rx_watch','statue','stillLife','rx_screws'};
% set_names = {'rx_motor','rx_screws'};
% set_names = {'scene1'};

for scale = 2:3
    disp(['----- Scale is ',num2str(scale)]);
    for ds = 1:2
        if ds == 1
            ds_method = 'bicubic';
        elseif ds == 2
            ds_method = 'gaussian';
        end
        disp(['------- Ds method is ',ds_method]);

        for set_ind = 1:length(set_names)
            %% settings
            set_name = set_names{set_ind};
            
            if set_ind <= 10 % for the last two, the dim is 7
                lf_dim = 9;
            else
                lf_dim = 7;
            end
            disp(['---------- Processing ',set_name,' with dim ',num2str(lf_dim)]);

            folder = ['/home/joechan/DataSet/HCI-dataset/',set_name,'.mat'];
            save_folder = ['/media/joechan/Standby/LFVDSR/Data_LFVDSR/train/HCI/train_LFVDSR_',set_name,'_scale',num2str(scale),'_',ds_method,'.mat'];

            size_input = 41;
            size_label = 41;

            stride = 30;

            %% initialization
            data = zeros(size_input, size_input, 4, 40000);
            label = zeros(size_label, size_label, 4, 40000);
            padding = abs(size_input - size_label)/2;
            count = 0;

            %% generate data

            load(folder);
%             lf_data = im2double(lf_data);
            seed_view = squeeze(lf_data(1,1,:,:,:));
            seed_view = modcrop(seed_view,scale);
            [H,W,~] = size(seed_view);

            lf_input = zeros([lf_dim,lf_dim,H,W],'uint8');
            lf_label = zeros([lf_dim,lf_dim,H,W],'uint8');

            for u = 1:lf_dim
                for v = 1:lf_dim
                    % get the label images after gray scale extraction and modcrop
                    sub_img = squeeze(lf_data(u,v,:,:,:));
                    sub_img = rgb2ycbcr(sub_img);
                    sub_img = squeeze(sub_img(:,:,1));
                    sub_img = modcrop(sub_img,scale);
                    lf_label(u,v,:,:) = sub_img;

                    % get the input images 
                    if strcmpi(ds_method,'bicubic')
                        sub_img_input = imresize(imresize(sub_img,1/scale,'bicubic'),scale,'bicubic');
                    elseif strcmpi(ds_method,'gaussian')
                        sub_img_input = imresize(downsample_gauss(sub_img,scale),scale,'bicubic');
                    else
                        disp('Wrong downsampling kernel name!');
                    end
                    lf_input(u,v,:,:) = sub_img_input;
                end
            end
            lf_input = im2double(lf_input);
            lf_label = im2double(lf_label);

            % now store them as patches
            for i = 1:lf_dim-1
                for j=1:lf_dim-1 
                    label1 = squeeze(lf_label(i,j,:,:));
                    label2 = squeeze(lf_label(i,j+1,:,:));
                    label3 = squeeze(lf_label(i+1,j,:,:));
                    label4 = squeeze(lf_label(i+1,j+1,:,:));

                    input1 = squeeze(lf_input(i,j,:,:));
                    input2 = squeeze(lf_input(i,j+1,:,:));
                    input3 = squeeze(lf_input(i+1,j,:,:));
                    input4 = squeeze(lf_input(i+1,j+1,:,:));

                    for x = 1:stride:(H-size_input+1)
                        for y = 1:stride:(W-size_input+1)
                            subim_label1 = label1(x : x+size_input-1, y : y+size_input-1);
                            subim_label2 = label2(x : x+size_input-1, y : y+size_input-1);
                            subim_label3 = label3(x : x+size_input-1, y : y+size_input-1);
                            subim_label4 = label4(x : x+size_input-1, y : y+size_input-1);

                            subim_input1 = input1(x+padding : x+padding+size_label-1, y+padding : y+padding+size_label-1);
                            subim_input2 = input2(x+padding : x+padding+size_label-1, y+padding : y+padding+size_label-1);
                            subim_input3 = input3(x+padding : x+padding+size_label-1, y+padding : y+padding+size_label-1);
                            subim_input4 = input4(x+padding : x+padding+size_label-1, y+padding : y+padding+size_label-1);

                            count = count + 1;
                            label(:,:,1,count) = subim_label1;
                            label(:,:,2,count) = subim_label2;
                            label(:,:,3,count) = subim_label3;
                            label(:,:,4,count) = subim_label4;

                            data(:,:,1,count) = subim_input1;
                            data(:,:,2,count) = subim_input2;
                            data(:,:,3,count) = subim_input3;
                            data(:,:,4,count) = subim_input4;
                        end
                    end
                end
            end
            

            data_1 = data(:,:,:,1:floor(count/2));
            data_2 = data(:,:,:,floor(count/2)+1:count);
            label_1 = label(:,:,:,1:floor(count/2));
            label_2 = label(:,:,:,floor(count/2)+1:count);

            save(save_folder,'data_1','data_2','label_1','label_2');
        end
    end
end