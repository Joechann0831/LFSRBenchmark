function SR_LF = bm_pca_rr_parent(LR_LF,mf)
% The bm_pca_rr_parent is the parent function that is used to configure the
% pca_rr light field super-resolution. This implementation is based on the
% work published in
%
%  R.A. Farrugia, C. Galea, C. Guillemot, "Super Resolution of Light Field 
%  Images using Linear Subspace Projection of Patch-Volumes," in IEEE 
%  Journal on Selected Topics in Signal Processing, vol. 11, no. 7, 
%  pp. 1058-1071, Oct. 2017
% 
% It takes the low-resolution light field LR_LF and the scale by which the
% light field will be restored mf as inputs and returns the restored light
% field. In summary the algorithm will extract patch volumes and use the
% pca basis and ridge regression matrices which are provided as a model to
% estimate the corresponding high resolution patch volume. More information
% is provided in the paper.
%
% Input: LR_LF: This corresponds to the low-resolution lightfield
%        param: structure containing different parameters, dictionaries and pca 
%               basis functions which were pre-computed
%
% Reuben Farrugia
% version 1 - 23/8/2016
% version 2 - 26/1/2018
%


addpath('./MATLAB/superresolution/pca_rr/');

%--------------------------------------------------------------------------
% Parameter configuration
%--------------------------------------------------------------------------
param.patch_size = 10;   % Patch size
param.Nl         = 500;  % Dimension of the low-resolution patch volume
param.Nh         = 500;  % Dimension of the high-resolution patch volume
param.Npts       = 2000; % Number of points from each light field
param.window_size= 8;    % Search window size 0 for pca_rr
param.lambda     = 1E-6; % Ridge Regression regularization paramete
param.overlap    = floor(param.patch_size/3);
%--------------------------------------------------------------------------
% Derive the filename where the pca basis will be stored
pca_foldername = './DATA/superresolution/pca_rr/';
pca_filename = sprintf('%sbm_pca_basis_x%d.mat',pca_foldername, mf);
if ~exist(pca_filename,'file')
    if ~exist(pca_foldername,'dir')
        % Create the foldername where to store the model
        mkdir(pca_foldername);
    end
    
    % Get the url where the file is stored
    if mf == 2
        %url = 'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=1SP0kGzN7kWRh-ZM2lAI_LXRQB51m2VLf&export=download';
        url = 'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=1LOhy-begjsXFad3sReQLpTRvxEEdDjb-&export=download';
    elseif mf == 3
        %url = 'https://drive.google.com/a/um.edu.mt/uc?authuser=1&id=1hQH-uIbiDxOIHAWMRRtuErmdM-XWeZTl&export=download';
        url = 'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=1qyo42A5zVcUSLo5MY6AEfmWzDltoA7UT&export=download';
    elseif mf == 4
        %url = 'https://drive.google.com/a/um.edu.mt/uc?authuser=1&id=1uBu01qHx_4quBhUh9OKjVOnrFfuJ0EqW&export=download';
        url = 'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=1smQmhKFZFnekrVp9IkDUzJFVaqaS7IYp&export=download';
    end
    % Download the required pca basis
    urlwrite(url,pca_filename);
end
% Load the PCA basis
data = load(pca_filename);

% Put the pca basis in the parameter
param.pca_basis = data.pca_basis;
param.phi       = data.phi;

% Clear the data
clearvars('data');

% Convert the low-resolution light field to ycbcr components
[LR_LF_y, LR_LF_cb,LR_LF_cr] = LF_rgb2ycbcr(LR_LF);

% Compute super-resolution using PCA+RR (search window was set to zero)
SR_LF_y = bm_pca_rr_LF(LR_LF_y, param);

% Convert the restored light field in the rgb color channel
SR_LF = LF_ycbcr2rgb(SR_LF_y,LR_LF_cb,LR_LF_cr);

