function pca_rr_training(method, mf)
% This script will be used to train the PCA+RR or BM+PCA+RR. The method can
% be either pca_rr or bm_pca_rr and considers magnification factors 2,3,4.
% The model will be restored in ../DATA/superresolution/pca_rr/
%

clc; close all;

addpath('../MATLAB/light_field/');
addpath('../MATLAB/general/');
addpath('../MATLAB/superresolution/pca_rr/');

param.Npts         = 1500;  % number of coordinate points
param.mf           = mf;     % magnification factor
param.patch_size   = 10;
param.Nl           = 500;
param.Nh           = 500;
lambda             = 1E-6;

if strcmp(method,'pca_rr')
    param.window_size = 0;
else
    param.window_size = 8;
end

% Load the light fields to be considered for training
[lf_names_test, ~] = read_configuration('../config\superresolution.cfg');

% Determine the list of light fields to be considered for training
[lf_names_train,datasets_train] = get_complement_list(lf_names_test,'../../LF-DATASET/');

% Determine the number of light fields to use for training
N = size(lf_names_train,2);

%--------------------------------------------------------------------------
% Choose random patch volumes
%--------------------------------------------------------------------------
% Initialize the dictionaries
Hd = zeros(8100,N*param.Npts,'uint8');
Ld = zeros(8100,N*param.Npts,'uint8');

for i = 1:N
    fprintf('Light Field %s %d out of %d\n',lf_names_train{i},i,N);
    
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

    % Generate the required number of patch coordinates
    patch_idx = get_patch_idx(size(HR_LF,1),size(HR_LF,2),param.patch_size,param.Npts);
    
    if param.window_size == 0
        % Get the high resolution dictionary with no disparity comp.
        Dlf_H = get_dictionary(HR_LF,patch_idx, param.patch_size);
        % Get the low resolution dictionary with no disparity comp.
        Dlf_L = get_dictionary(LR_LF,patch_idx, param.patch_size);
    else
        % Get the coupled dictionaries aligned using block matching
        [Dlf_H, Dlf_L] = get_dictionary_aligned(LR_LF, HR_LF,patch_idx,param.patch_size, param.window_size); 
    end
    
    % Determine the index values where the patch volumes will be stored
    strt_idx = (i-1)*param.Npts+1;
    end_idx  = i*param.Npts;

    % Put the retrieved patch volumes in the dictionary
    Hd(:,strt_idx:end_idx) = Dlf_H;
    Ld(:,strt_idx:end_idx) = Dlf_L;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Learn PCA basis
%--------------------------------------------------------------------------
pca_basis = learn_pca_basis(Ld, Hd, param);

% Put the mean vector for the low and high resolution patch volumes in the
% output structure
Ml = pca_basis.Ml; Mh = pca_basis.Mh;
El = pca_basis.El; Eh = pca_basis.Eh;

% Convert points to single precision
Ld = single(Ld); Hd = single(Hd);

% Project the low-resolution dictionary on the sub-space
Ls = El' * (Ld - repmat(Ml,[1,size(Ld,2)]));

% Project the high-resolution dictionary on the sub-space
Hs = Eh' * (Hd - repmat(Mh,[1,size(Hd,2)]));

% Compute the upscaling function using all elements in the dictionary. This 
% projection matrix will be used to project a point from the low-quality onto 
% the high quality sub-spaces
phi = Hs * Ls' /(Ls*Ls' + lambda*eye(size(Ls,1)));

% Derive the upscaling function psi
%psi = Eh*phi*El';
%--------------------------------------------------------------------------

if strcmp(method,'pca_rr')
    % Derive the output filename
    out_filename = sprintf('../DATA/superresolution/pca_rr/pca_basis_x%d.mat',mf);
elseif strcmp(method,'bm_pca_rr')
    % Derive the output filename
    out_filename = sprintf('../DATA/superresolution/pca_rr/bm_pca_basis_x%d.mat',mf);   
end

% Put pca-basis into a structure
pca_basis.Eh = Eh; pca_basis.El = El; pca_basis.Mh = Mh; pca_basis.Ml = Ml;

% Save the trained model
save(out_filename,'phi','pca_basis');