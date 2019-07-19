function [Dhf, Dlf] = get_dictionary_aligned(LR_LF, HR_LF,patch_idx,patch_size, window_size)

% Choose the anchor sub-aparture image
anchor_sai_idx = median(1:size(LR_LF,3));

% Initialize the low and high resolution dictionaries
Dlf = zeros(patch_size*patch_size*size(LR_LF,3)*size(LR_LF,4), size(patch_idx,1));
Dhf = zeros(patch_size*patch_size*size(LR_LF,3)*size(LR_LF,4), size(patch_idx,1));

for i = 1:size(patch_idx,1)
    %fprintf('%d out of %d patches included in dictionary\n',i,size(patch_idx,1));
    %fflush(stdout);
    % Get the region from which the patch will be extracted
    n_patch_strt = patch_idx(i,1); 
    m_patch_strt = patch_idx(i,2);
    n_patch_end  = n_patch_strt + patch_size-1; 
    m_patch_end  = m_patch_strt + patch_size-1;
    
    % Extract the current patch to consider
    x_patch = get_anchor_patch(LR_LF,n_patch_strt,n_patch_end, m_patch_strt,m_patch_end, ...
                               anchor_sai_idx);
        
    % Derive the search window
    search_window = get_search_window(n_patch_strt, n_patch_end, m_patch_strt, m_patch_end,size(LR_LF),window_size);
    
    % Derive the low-resolution aligned patch volume
    [LR_patch_volume,dv] = get_LR_aligned_patch_volume(LR_LF,search_window,x_patch,patch_size,anchor_sai_idx);

    % Initialize the HR patch volume
    HR_patch_volume = zeros(size(LR_patch_volume));
    for j = 1:size(HR_LF,3)
       if j == anchor_sai_idx
            idx_m = m_patch_strt:min([m_patch_end,size(HR_LF,2)]);
            idx_n = n_patch_strt:min([n_patch_end,size(HR_LF,1)]);
        else
            idx_m = search_window{1}.m_strt + dv(j,2) - 1:min([search_window{1}.m_strt+patch_size-1+dv(j,2)-1,size(HR_LF,2)]);
            idx_n = search_window{1}.n_strt + dv(j,1) - 1:min([search_window{1}.n_strt+patch_size-1+dv(j,1)-1,size(HR_LF,1)]);
        end
        
        % Put the result into HR_patch
        HR_patch_volume(:,:,j) = HR_LF(idx_n,idx_m,j);
    end
    % Extract the vector representing the low-resolution aligned LF
    x = reshape(LR_patch_volume, [patch_size*patch_size*size(LR_LF,3),1]);
    % Extract the vector representing the high-resolution aligned LF
    y = reshape(HR_patch_volume, [patch_size*patch_size*size(LR_LF,3),1]);

    % Put the current volume patch
    Dlf(:,i) = x; Dhf(:,i) = y;
end



