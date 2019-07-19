function idx = get_patch_idx(Hlf,Wlf,patch_size,Npts)
% This function returns the patch index of the chosen patches. The patches
% are randomly selected
idx = [];

k = 1;
while 1
    % Generate a random patch coordinate
    y_strt = randi([1 Hlf],1,1); x_strt = randi([1 Wlf],1,1);
    y_end  = y_strt + patch_size - 1; x_end = x_strt + patch_size;
    
    if y_end <= Hlf && x_end <= Wlf % This is a valid patch
        % Put this in index
        idx = [idx; y_strt, x_strt];
        % Increment the size
        k = k + 1;
    end
    
    if k > Npts*2
        break;
    end
end

idx = unique(idx,'rows');

idx = idx(1:Npts,:);
