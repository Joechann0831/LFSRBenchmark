function LF_BLINDdenoising_analysis(dn_method,sig,out_flag)
% The LF_denoising_analysis script is a command line script that
% allows to analyse the performance of a number of light field
% denoising algorithms.
% 
% Syntax: LF_denoising_analysis('bm3d',10)
%
% Input: dn_method: which specifies the super-resolution algorithm to be 
%                   simulated. The following is a list of sr_methods that 
%                   are supported here:
%                   - 'bm3d' - classical BM3D denoising algorithm of each
%                                 sub-aperture image independently [1]
%                   - 'vbm4d' - extension of the bm3d algorithm where they
%                   consider denoising of 4D video [2]
%
%        sig: standard deviation of the noise
%
%        out_flag: This is a boolean value which specifies if the results
%        attined in the simulation will be stored or not. By default it is
%        set to false
%
% [1] K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, “Image Denoising
% by Sparse 3-D Transform-Domain Collaborative Filtering,” IEEE Trans.
% on Image Processing, vol. 16, no. 8, pp. 2080–2095, 2007
%
% [2] M. Maggioni, G. Boracchi, A. Foi, and K. Egiazarian, “Video denoising,
% deblocking, and enhancement through separable 4-D nonlocal spatiotemporal
% transforms,” IEEE Trans. on Image Processing, vol. 21, no. 9, pp. 
% 3952–3966, 2012
% 
% [3] Z. Li, H. Baker, and R. Bajcsy, “Joint image denoising using light-field
% data,” Proc. ICMEW, 2013. 
%

clc; close all;
addpath('MATLAB/');
addpath('MATLAB/light_field/');
addpath('MATLAB/denoising/');

warning off; rng(0);

%--------------------------------------------------------------------------
% SIMULATION CONFIGURATION
%--------------------------------------------------------------------------

if nargin == 2
    out_flag = false;
end

% Load the light fields to be processed in this simulation
[lf_names, datasets] = read_configuration('config\BLINDdenoising.cfg');

% Determine the number of light fields to consider
N = size(lf_names,2);

fprintf('--------------------------------------------------------------\n');
if strcmp(dn_method,'noisy')
    fprintf('Evaluating the performance of Noisy\n');
    out_filename = sprintf('RESULTS/BLINDdenoising/sig%d/noisy.csv',sig);
    out_img_foldername = sprintf('RESULTS/BLINDdenoising/sig%d/centre_view/noisy/',sig);
    out_LF_foldername  = sprintf('RESULTS/BLINDdenoising/sig%d/LF/noisy/',sig);
elseif strcmp(dn_method,'bm3d')
    fprintf('Evaluating the performance of BM3D\n');
    out_filename = sprintf('RESULTS/BLINDdenoising/sig%d/bm3d.csv',sig);
    out_img_foldername = sprintf('RESULTS/BLINDdenoising/sig%d/centre_view/bm3d/',sig);
    out_LF_foldername  = sprintf('RESULTS/BLINDdenoising/sig%d/LF/bm3d/',sig);
elseif strcmp(dn_method,'vbm4d')
    fprintf('Evaluating the performance of VBM4D\n');
    out_filename = sprintf('RESULTS/BLINDdenoising/sig%d/vbm4d.csv',sig);
    out_img_foldername = sprintf('RESULTS/BLINDdenoising/sig%d/centre_view/vbm4d/',sig);
    out_LF_foldername  = sprintf('RESULTS/BLINDdenoising/sig%d/LF/vbm4d/',sig);
elseif strcmp(dn_method,'bm3d-epi')
    fprintf('Evaluating the performance of BM3D-EPI\n');
    out_filename = sprintf('RESULTS/BLINDdenoising/sig%d/bm3d-epi.csv',sig);
    out_img_foldername = sprintf('RESULTS/BLINDdenoising/sig%d/centre_view/bm3d-epi/',sig);
    out_LF_foldername  = sprintf('RESULTS/BLINDdenoising/sig%d/LF/bm3d-epi/',sig);
elseif strcmp(dn_method,'lfbm5d')
    fprintf('Evaluating the performance of LFBM5D\n');
    out_filename = sprintf('RESULTS/BLINDdenoising/sig%d/lfbm5d.csv',sig);
    out_img_foldername = sprintf('RESULTS/BLINDdenoising/sig%d/centre_view/lfbm5d/',sig);
    out_LF_foldername  = sprintf('RESULTS/BLINDdenoising/sig%d/LF/lfbm5d/',sig);
elseif strcmp(dn_method,'b-bm3d')
    fprintf('Evaluating the performance of B-BM3D\n');
    out_filename = sprintf('RESULTS/BLINDdenoising/sig%d/b-bm3d.csv',sig);
    out_img_foldername = sprintf('RESULTS/BLINDdenoising/sig%d/centre_view/b-bm3d/',sig);
    out_LF_foldername  = sprintf('RESULTS/BLINDdenoising/sig%d/LF/b-bm3d/',sig);
end
fprintf('--------------------------------------------------------------\n');

% Make sure that the folder containing the images is available
if ~exist(out_img_foldername,'dir')
    mkdir(out_img_foldername);
end

% Make sure that the folder containing the light fields is available
if ~exist(out_LF_foldername,'dir')
    mkdir(out_LF_foldername);
end

% Open the file where the results will be stored
if out_flag
    out_filename_st = out_filename;
    while exist(out_filename,'file')
        out_filename = out_filename_st;
        if exist(out_filename,'file')
            out_filename2 = [out_filename(1:end-4),sprintf('%d',randi(1000000)),'.csv'];
        end
    	out_filename = out_filename2;
    end
    
    fid = fopen(out_filename,'w');
end

for n = 1:N
    % Get the dataset
    dataset = datasets{n};
    % Get the light field name
    lf_name = lf_names{n};
    
    out_img_filename = [out_img_foldername,lf_name,'.bmp'];
    
    out_LF_filename  = [out_LF_foldername,lf_name,'.mat'];
    
    % Derive the folder containing the light field
    dataset_foldername = sprintf('../LF-DATASET/%s/',dataset);
    
    % Load the high resolution light field
    if strcmp(dataset,'HCI')
        % Load the high resolution light field
        HQ_LF = load_hci_lf(dataset_foldername, lf_name);
    elseif strcmp(dataset,'STANFORD')
        HQ_LF = load_stanford_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'EPFL')
        HQ_LF = load_epfl_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'INRIA')
        HQ_LF = load_inria_lf(dataset_foldername,lf_name);
    end
    % Permute the dimensions
    HQ_LF = permute(HQ_LF,[3,4,5,1,2]);

    % Degrade the input lighfield using AWGN
    LQ_LF = lf_awgn(HQ_LF,sig);
    
    if strcmp(dn_method,'noisy')
        RQ_LF = LQ_LF;
    elseif strcmp(dn_method,'bm3d')
        % Use bm3d denoising
        RQ_LF = LF_bm3d_denoising(LQ_LF);
    elseif strcmp(dn_method,'vbm4d')
        % Use vbm4d denoising
        RQ_LF = LF_vbm4d_denoising(LQ_LF);
    elseif strcmp(dn_method,'bm3d-epi')
        % Use the BM3d-EPI method
        RQ_LF = LF_bm3d_epi_denoising(LQ_LF);
    elseif strcmp(dn_method,'lfbm5d')
        % Use the LFBM5D method
        RQ_LF = LFBM5D_denoising(LQ_LF);
    elseif strcmp(dn_method,'b-bm3d')
        RQ_LF = LF_B_bm3d(LQ_LF);
    end
            
    % Extract the centre view
    Ic = RQ_LF(:,:,:,5,5);
    
    if out_flag
        % Write the center view
        imwrite(Ic,out_img_filename);
    end
    
    % Compute the psnr evaluation
    [psnr_val,ssim_val] = quality_analysis(HQ_LF, RQ_LF);
    
    fprintf('%s, %s, %0.4f, %0.4f\n',lf_name,dataset,psnr_val, ssim_val);
    if out_flag
        fprintf(fid,'%s, %s, %0.4f, %0.4f\n',lf_name,dataset,psnr_val, ssim_val);
    
        % Save the restored light field
        save(out_LF_filename,'RQ_LF');
    end
end
fclose('all');

