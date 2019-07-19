function [x_patch_,mv,error] = block_match(window_i, x_patch)
patch_size = size(x_patch,1);
window_i = double(window_i);
x_patch  = double(x_patch);
error = inf;
for i = 1:size(window_i,1)
    for j = 1:size(window_i,2)
        if i+patch_size-1 <= size(window_i,1) && j+patch_size-1 <= size(window_i,2)
            x_patch_ij = window_i(i:i+patch_size-1, j:j+patch_size-1);    
            eij = sum(sum(abs(x_patch_ij - x_patch)));
            if eij < error
                x_patch_ = x_patch_ij;
                mv = [i,j];
                error = eij;
            end
        end
    end
end
x_patch_ = uint8(x_patch_);
% Normalize the error per pixel
error = error/patch_size^2;