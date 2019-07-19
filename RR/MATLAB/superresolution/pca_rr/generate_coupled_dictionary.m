function [D_hrlf,D_lrlf] = generate_coupled_dictionary(param,list,Npts)

D_hrlf = []; D_lrlf = [];
for i = 1:size(list,1)
    % This image is used for training
    lf_name = list{i}.filename;
    dataset = list{i}.dataset;
    
    dataset_foldername = sprintf('DATASET/%s/',dataset);
    if strcmp(dataset,'INRIA')
        HR_LF = load_inria_lf(dataset_foldername,lf_name(1:end-4));
    elseif strcmp(dataset,'STANFORD')
        HR_LF = load_stanford_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'HCI')
        HR_LF = load_hci_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'EPFL')
        HR_LF = load_epfl_lf(dataset_foldername,lf_name(1:end-4));
    end

    % Generate low quality light field
    LR_LF = generate_low_quality_lightfield(HR_LF,param);

    % Decompose the light field into luminance and chrominance color channels
    [HR_LF, ~, ~] = LF_rgb2ycbcr(HR_LF);
    [LR_LF, ~, ~] = LF_rgb2ycbcr(LR_LF);
      
    % Generate the required number of patch coordinates
    %fprintf('Get the indices of the random patches centered on the center view.\n');
    patch_idx = get_patch_idx(size(HR_LF,1),size(HR_LF,2),param.pca_rr.patch_size,Npts);
    
    if param.pca_rr.window_size == 0
        % Get the high resolution dictionary with no disparity comp.
        Dlf_H = get_dictionary(HR_LF,patch_idx, param.pca_rr.patch_size);
        % Get the low resolution dictionary with no disparity comp.
        Dlf_L = get_dictionary(LR_LF,patch_idx, param.pca_rr.patch_size);
    else
        % Get the coupled dictionaries aligned using block matching
        [Dlf_H, Dlf_L] = get_dictionary_aligned(LR_LF, HR_LF,patch_idx,param.pca_rr.patch_size, param.pca_rr.window_size); 
    end
    % Save the dictionary of the high-resolution lightfield
    D_hrlf = [D_hrlf,uint8(Dlf_H)];
    D_lrlf = [D_lrlf,uint8(Dlf_L)];
    
    fprintf('%d out of %d Lfs\n',i,size(list,1));
end
