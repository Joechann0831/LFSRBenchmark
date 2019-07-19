function [patch_volume,dv] = get_LR_aligned_patch_volume(LF,search_window,x_patch,patch_size,anchor_sai_idx)

% Initialize the patch volume to negative ones
patch_volume = zeros(patch_size,patch_size,size(LF,3));

% Derive the aligned patch volume
dv = zeros(size(LF,3),2);

for j = 1:size(LF,3)
    if j ~= anchor_sai_idx
        % Get the window where to search
        window_i = get_search_window_sub_ap_img(LF,search_window{1},j,patch_size);
        % Use block-match to find the closest patch to x_patch
        [x_patch_,dv(j,:)] = block_match(window_i, x_patch);
        % Put the aligned patch into the volume patch
        patch_volume(:,:,j) = x_patch_;        
    else
        % Put the anchor patch into volume patch
        patch_volume(:,:,j) = x_patch;
    end
end



