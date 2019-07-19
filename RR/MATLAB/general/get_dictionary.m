function  Dlf = get_dictionary(L,idx, patch_size)
% Initialize the dictionary
Dlf = zeros(patch_size*patch_size*size(L,3)*size(L,4), size(idx,1));
for i = 1:size(idx,1)
    % Get the start and end coordinates for the patch
    y_strt = idx(i,1); x_strt = idx(i,2);
    y_end = y_strt+patch_size-1; x_end = x_strt+patch_size-1;
    
    % extract the volume from the lightfield
    x = reshape(L(y_strt:y_end,x_strt:x_end,:,:), [patch_size*patch_size*size(L,3)*size(L,4),1]);
    
    % Put the current volume patch
    Dlf(:,i) = x;
end

