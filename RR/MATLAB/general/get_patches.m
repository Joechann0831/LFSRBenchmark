function  Dlf = get_patches(L,idx, patch_size)
% Initialize the dictionary
Dlf = zeros(patch_size,patch_size, 1, size(idx,1),'uint8');
for i = 1:size(idx,1)
    % Get the start and end coordinates for the patch
    y_strt = idx(i,1); x_strt = idx(i,2);
    y_end = y_strt+patch_size-1; x_end = x_strt+patch_size-1;
    
    % extract the volume from the lightfield
    x = L(y_strt:y_end,x_strt:x_end);
    
    % Put the current volume patch
    Dlf(:,:,1,i) = x;
end

