function get_random_patches_lfsrcnn(mf)
% This function is used to train the light field SRCNN method using random
% patches from a set of training data. This function will store the imdb
% structure that will be used by matconvnet to train the LFSRCNN network.
%

clc; close all;

addpath('../MATLAB/light_field/');
addpath('../MATLAB/general/');

param.Npts         = 150;   % number of coordinate points from each SAI
param.mf           = mf;    % magnification factor
param.patch_size   = 32;    % patch size

% Load the light fields to be considered for training
[lf_names_test, ~] = read_configuration('../config\superresolution.cfg');

% Determine the list of light fields to be considered for training
[lf_names_train,datasets_train] = get_complement_list(lf_names_test,'../../LF-DATASET/');

% Determine the number of
N = size(lf_names_train,2);

% Initialize the imdb data where the patch-volumes will be stored
imdb.images.data = zeros(param.patch_size,param.patch_size, 1 ,N*param.Npts*81,'uint8');
imdb.images.labels = zeros(param.patch_size,param.patch_size,1,N*param.Npts*81,'uint8');
imdb.images.set = 2*ones(1,N*param.Npts*81);

% Determine the number of
M = size(imdb.images.set,2);
Ntrain = round(0.8*M);
idx = randperm(M,Ntrain);

imdb.images.set(idx) = 1;             % set 1 train data
imdb.images.id(1,:) = 1:size(imdb.images.data,4);   % give unique ids 

j = 1;
for i = 1:N
    fprintf('Light Field %s %d out of %d\n',lf_names_train{i},i,N);
    %----------------------------------------------------------------------
    % Load data
    %----------------------------------------------------------------------
    dataset = datasets_train{i};
    
    fprintf('  Load high-resolution light-field\n');
    dataset_foldername = sprintf('../../LF-DATASET/%s/',dataset);
    lf_name = lf_names_train{i};
    
    % Load the high-resolution light-field
    if strcmp(dataset,'HCI')
        % Load the high resolution light field
        HR_LF = load_hci_lf(dataset_foldername, lf_name);
    elseif strcmp(dataset,'STANFORD')
        HR_LF = load_stanford_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'EPFL')
        HR_LF = load_epfl_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'INRIA')
        HR_LF = load_inria_lf(dataset_foldername,lf_name);
    end
    fprintf('  Generate low-resolution light-field\n');
    
    HR_LF = permute(HR_LF,[3,4,5,1,2]);
    
    % Generate the low-resolution light field
    LR_LF = lf_downsample(HR_LF,mf);
    
    % Decompose the light field into luminance and chrominance color channels
    [HR_LF, ~, ~] = LF_rgb2ycbcr(HR_LF);
    [LR_LF, ~, ~] = LF_rgb2ycbcr(LR_LF);

    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Choose random patches
    %----------------------------------------------------------------------
    for k = 1:size(LR_LF,3)
        % Generate the required number of patch coordinates
        patch_idx = get_patch_idx(size(LR_LF,1),size(LR_LF,2),param.patch_size,param.Npts);
        
        % Determine the index values where the patch volumes will be stored
        strt_idx = (j-1)*param.Npts+1;
        if k == 1
            strt_idx1 = strt_idx;
        end
        end_idx  = j*param.Npts;
        
        % Get the high resolution dictionary with no disparity comp.
        imdb.images.labels(:,:,:,strt_idx:end_idx) = uint8(get_patches(HR_LF(:,:,k),patch_idx, param.patch_size));
        % Get the low resolution dictionary with no disparity comp.
        imdb.images.data(:,:,:,strt_idx:end_idx) = uint8(get_patches(LR_LF(:,:,k),patch_idx, param.patch_size));

        j = j + 1;
    end
    fprintf(1,'Strt %d End %d\n',strt_idx1,end_idx);
    %----------------------------------------------------------------------
end
imdb_foldername = '../DATA/superresolution/lfsrcnn/';

if ~exist(imdb_foldername,'dir')
    mkdir(imdb_foldername);
end

% Define the imdb structure where the data will be stored
imdb_filename = sprintf('../DATA/superresolution/lfsrcnn/imdb_x%d.mat',mf);

% Save the imdb data
fprintf('Save the imdb structure...\n');
save(imdb_filename,'imdb','-v7.3');
