%%% This script is used to write the training data to h5 files

clc;
clear;

%% settings
whole_set = {'buddha','buddha2','coneHead','medieval','monasRoom','papillon','rx_elephant','rx_watch','statue','stillLife','rx_screws'};
test_set = {'buddha','buddha2','monasRoom'};
test_mode = 1;
scale = 2;
ds_method = 'bicubic';

input_size = 41;
label_size = 41;


%% load the data and write them to H5 files
% train_set = setdiff(whole_set,test_set);
train_set = test_set;

folder_prefix = '/media/joechan/Standby/LFVDSR/Data_LFVDSR/train/HCI/train_LFVDSR_';
folder_suffix = ['_scale',num2str(scale),'_',ds_method,'.mat'];

diary('./train_seq_HCI_scale2_bicubic.txt');
diary on;

for ind = 1:length(train_set)
    data = zeros(input_size,input_size,4,50000);
    label = zeros(input_size,input_size,4,50000);
    count = 0;
    
    MatData = load([folder_prefix,train_set{ind},folder_suffix]);
    [~,~,~,num1] = size(MatData.data_1);
    [~,~,~,num2] = size(MatData.data_2);
    disp(['patch number: ',num2str(num1)]);
    disp(['patch number: ',num2str(num2)]);
    
    data(:,:,:,count+1:count+num1) = MatData.data_1;
    label(:,:,:,count+1:count+num1) = MatData.label_1;
    count = count + num1;
    data(:,:,:,count+1:count+num2) = MatData.data_2;
    label(:,:,:,count+1:count+num2) = MatData.label_2;
    count = count + num2;
    
    order = randperm(count);
    data = data(:,:,:,order);
    label = label(:,:,:,order);
    
    savepath = ['/media/joechan/Joechan/LFVDSR/Training/HCI/scale',num2str(scale),'/',ds_method,'/train_data_',train_set{ind},'.h5'];
    
    chunksz = 64;
    created_flag = false;
    totalct = 0;
        
    for batchno = 1:floor(count/chunksz)
        last_read=uint64((batchno-1)*chunksz);
        batchdata = data(:,:,:,last_read+1:last_read+chunksz);
        batchlabs = label(:,:,:,last_read+1:last_read+chunksz);%batchlabs = double(label(:,:,1,last_read+1:last_read+chunksz))/255.0;
        startloc = struct('dat',[1,1,1,totalct+1], 'lab', [1,1,1,totalct+1]);
        curr_dat_sz = store2hdf5(savepath, batchdata, batchlabs, ~created_flag, startloc, chunksz);
        created_flag = true;
        totalct = curr_dat_sz(end);
    end
    h5disp(savepath);
end

diary off;

