function [y, mask] = stitching(x_prime,mask,x)
% Clip the recovered patch to ensure it is within range
x = min(max(x,0),255);
% concatenate to mask to make it same size as x
msk = -1*ones(size(x));
msk(1:size(x_prime,1),1:size(x_prime,2)) = mask;
mask = msk;
xx_prime = zeros(size(x));
xx_prime(1:size(x_prime,1),1:size(x_prime,2)) = x_prime;
y = -ones(size(x));
y(mask == 0) = x(mask == 0);
y(mask == 1) = (x(mask == 1) + xx_prime(mask == 1))/2;
y    = y(1:size(x_prime,1),1:size(x_prime,2));