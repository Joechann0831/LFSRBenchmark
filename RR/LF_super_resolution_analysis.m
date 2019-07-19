function LF_super_resolution_analysis(sr_method,mf,out_flag)
% The LF_super_resolution_analysis script is a command line script that
% allows to analyse the performance of a number of light field
% super-resolution algorithms.
% 
% Syntax: LF_super_resolution_analysis('bicubic',3)
%
% Input: sr_method: which specifies the super-resolution algorithm to be 
%                   simulated. The following is a list of sr_methods that 
%                   are supported here:
%                   - 'bicubic' - classical bicubic interpolation of each
%                                 sub-aperture image independently
%                   - 'lf_srcnn' - this method applies a CNN to restore
%                                  each sub-ape[rture image separately from 
%                                  the others [1],[2]
%                   - 'vdsr' - this method applies the VDSR single image
%                             super reslution algorithm to restore each
%                             sub-aperture image separately from the
%                             others [3]. 
%                   - 'pca_rr' - this method applies the pca_rr published
%                                in [4] - super-resolves patch volumes
%
%                   - 'graph-lfsr' - uses the graph-based SR method in [5]
%
%                   - 'bm_pca_rr' - this method applies the pm_pca_rr
%                                   published in [4] - super-resolves 
%                                   patch volumes
%                   - 'pb-srcnn' - this method is the proposed method using
%                                  SRCNN to restore the principal basis
%                   - 'pb-vdsr'  - this method is the proposed method using
%                                  VDSR to restore the principal basis
%                   - 'pb-Lab402'- this method is proposed using Lab402 to
%                                  restore the principal basis
%
%        mf: numeric value that stands for the magnification factor that
%        the method has to super-resolve
%
%        out_flag: This is a boolean value which specifies if the results
%        attined in the simulation will be stored or not. By default it is
%        set to false
%
% [1] Y. Yoon, H. G. Jeon, D. Yoo, J. Y. Lee and I. S. Kweon, "Light-Field Image 
% Super-Resolution Using Convolutional Neural Network," in IEEE Signal 
% Processing Letters, vol. 24, no. 6, pp. 848-852, June 2017.
%                      
% [2] C. Dong, C. C. Loy, K. He and X. Tang, "Image Super-Resolution Using 
% Deep Convolutional Networks," in IEEE Transactions on Pattern Analysis 
% and Machine Intelligence, vol. 38, no. 2, pp. 295-307, Feb. 1 2016.
% 
% [3] J. Kim, J. K. Lee and K. M. Lee, "Accurate Image Super-Resolution Using 
% Very Deep Convolutional Networks," 2016 IEEE Conference on Computer Vision 
% and Pattern Recognition (CVPR), Las Vegas, NV, 2016, pp. 1646-1654.
% 
% [4] R.A. Farrugia, C. Galea, C. Guillemot, "Super Resolution of Light Field 
% Images using Linear Subspace Projection of Patch-Volumes," in IEEE 
% Journal on Selected Topics in Signal Processing, vol. 11, no. 7, 
% pp. 1058-1071, Oct. 2017
%
% [5] M. Rossi and P. Frossard, "Graph-based light field super-resolution," 
% 2017 IEEE 19th International Workshop on Multimedia Signal Processing 
% (MMSP), Luton, 2017, pp. 1-6
%
clc; close all;

addpath('MATLAB/');
addpath('MATLAB/light_field/');
addpath('MATLAB/superresolution/');

warning off;

%--------------------------------------------------------------------------
% SIMULATION CONFIGURATION
%--------------------------------------------------------------------------

if nargin == 2
    out_flag = false;
end
% Load the light fields to be processed in this simulation
[lf_names, datasets] = read_configuration('config\superresolution.cfg');

if strcmp(sr_method,'srcnn') || strcmp(sr_method,'lf_srcnn')
    addpath('MATLAB/superresolution/srcnn/');
elseif strcmp(sr_method,'vdsr')
    addpath('MATLAB/superresolution/vdsr/');
end
% Determine the number of light fields to consider
N = size(lf_names,2);

fprintf('--------------------------------------------------------------\n');
if (strcmp(sr_method,'srcnn') || (strcmp(sr_method,'lf_srcnn')))
    fprintf('Evaluating the performance of SRCNN\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/srcnn.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/srcnn/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/srcnn/',mf);
elseif strcmp(sr_method,'bicubic')
    fprintf('Evaluating the performance of Bicubic interpolation\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/bicubic.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/bicubic/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/bicubic/',mf);
elseif  strcmp(sr_method,'vdsr')
    fprintf('Evaluating the performance of VDSR\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/vdsr.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/vdsr/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/vdsr/',mf);
elseif strcmp(sr_method,'pca_rr')
    fprintf('Evaluating the performance of PCA+RR\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/pca_rr.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/pca_rr/',mf);
    out_LF_foldername = sprintf('RESULTS/superresolution/x%d/LF/pca_rr/',mf);
elseif strcmp(sr_method,'bm_pca_rr')
    fprintf('Evaluating the performance of BM+PCA+RR\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/bm_pca_rr.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/bm_pca_rr/',mf);
    out_LF_foldername = sprintf('RESULTS/superresolution/x%d/LF/bm_pca_rr/',mf);
elseif strcmp(sr_method,'pb-srcnn')
    fprintf('Evaluating the performance of PB-SRCNN\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/pb-srcnn.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/pb-srcnn/',mf);
    out_LF_foldername = sprintf('RESULTS/superresolution/x%d/LF/pb-srcnn/',mf);
elseif strcmp(sr_method,'pb-vdsr')
    fprintf('Evaluating the performance of PB-VDSR\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/pb-vdsr.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/pb-vdsr/',mf);
    out_LF_foldername = sprintf('RESULTS/superresolution/x%d/LF/pb-vdsr/',mf);
elseif strcmp(sr_method,'pb-lab402')
    % Initialize the gpu
    g=gpuDevice(1);
    reset(g); %GPU reset

    fprintf('Evaluating the performance of PB-Lab402\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/pb-lab402.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/pb-lab402/',mf);
    out_LF_foldername = sprintf('RESULTS/superresolution/x%d/LF/pb-lab402/',mf);
elseif strcmp(sr_method,'graph-lfsr')
    fprintf('Evaluating the performance of GRAPH-SR\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/graph-sr.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/graph-sr/',mf);
    out_LF_foldername = sprintf('RESULTS/superresolution/x%d/LF/graph-sr/',mf);
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
        HR_LF = load_hci_lf(dataset_foldername, lf_name);
    elseif strcmp(dataset,'STANFORD')
        HR_LF = load_stanford_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'EPFL')
        HR_LF = load_epfl_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'INRIA')
        HR_LF = load_inria_lf(dataset_foldername,lf_name);
    end
    
    % Permute the dimensions
    HR_LF = permute(HR_LF,[3,4,5,1,2]);
    
    if strcmp(sr_method,'pb-lab402') || strcmp(sr_method,'bicubic') ...
            || strcmp(sr_method,'graph-lfsr')
        % Generate the low-resolution light field
        LR_LF = lf_downsample(HR_LF,mf,1);
    else
        % Generate the low-resolution light field
        LR_LF = lf_downsample(HR_LF,mf);
    end        
    % The light field will be restored one sub-aperture image at a time
    % and will not exploit the light field structure
    if strcmp(sr_method,'bicubic')
        % Bicubic interpolation of each sub-aperture image separately
        SR_LF = lf_bicubic(LR_LF,size(HR_LF),mf);
    elseif strcmp(sr_method,'srcnn') || strcmp(sr_method,'lf_srcnn')
        % Super-resolution using LF-SRCNN[1],[2]
        SR_LF = lf_srcnn(LR_LF,mf);
    elseif strcmp(sr_method,'vdsr')
        % Super-resolution using VDSR on each sub-aperture image [3]
        SR_LF = lf_vdsr(LR_LF,mf);
    elseif strcmp(sr_method,'pca_rr')
        % Super-resolution using PCA+RR [4]
        SR_LF = pca_rr_parent(LR_LF,mf);
    elseif strcmp(sr_method,'bm_pca_rr')
        % Super-resolution using BM+PCA+RR
        SR_LF = bm_pca_rr_parent(LR_LF,mf);
    elseif strcmp(sr_method,'pb-srcnn')
        % Super-resolve using principal basis light field super resolution
        % using SRCNN
        SR_LF = principal_basis_SR(LR_LF,mf,lf_name,'srcnn');
    elseif strcmp(sr_method,'pb-vdsr')
        % Super-resolve using principal basis light field super resolution
        % using VDSR
        SR_LF = principal_basis_SR(LR_LF,mf,lf_name,'vdsr');
    elseif strcmp(sr_method,'pb-lab402')
        % Super-resolve using principal basis light field super-resolution
        % using lab402 method
        SR_LF = principal_basis_SR(LR_LF,mf,lf_name,'lab402',[size(HR_LF,1),size(HR_LF,2)]);
    elseif strcmp(sr_method,'graph-lfsr')
        % Super-resoluve using the graph based super-resolution method
        SR_LF = graph_based_SR(LR_LF,mf,size(HR_LF));
    end
            
    % Extract the centre view
    Ic = SR_LF(:,:,:,5,5);
    
    if out_flag
        % Write the center view
        imwrite(Ic,out_img_filename);
    end
    
    % Compute the psnr evaluation
    [psnr_val,ssim_val] = quality_analysis(HR_LF, SR_LF);
    
    fprintf('%s, %s, %0.4f, %0.4f\n',lf_name,dataset,psnr_val, ssim_val);
    if out_flag
        fprintf(fid,'%s, %s, %0.4f, %0.4f\n',lf_name,dataset,psnr_val, ssim_val);
    
        % Save the restored light field
        save(out_LF_filename,'SR_LF');
    end
end
fclose('all');

