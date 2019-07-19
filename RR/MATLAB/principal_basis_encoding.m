function [Ipb,param] = principal_basis_encoding(B_r, B_g, B_b,sz)

tol = 1E-5;
% Initialize a tensor to contain the principal basis
Ipb = zeros(sz(1),sz(2),3);

% Compute min-max normalization
[B_r, min_r, max_r] = minmax_normalization(B_r);
[B_g, min_g, max_g] = minmax_normalization(B_g);
[B_b, min_b, max_b] = minmax_normalization(B_b);

% Reshape the RGB components of the principal basis
B_r = reshape(B_r,[sz(1),sz(2)]);
B_g = reshape(B_g,[sz(1),sz(2)]);
B_b = reshape(B_b,[sz(1),sz(2)]);

% Reconstruct the principal basis in rgb
Ipb(:,:,1) = B_r;
Ipb(:,:,2) = B_g;
Ipb(:,:,3) = B_b;

% Ensure that we show the principal component (not its complement)
if min_r(1) < 0 && abs(min_r(1)) > tol
    Ipb(:,:,1) = 1 - Ipb(:,:,1);
end
if min_g(1) < 0 && abs(min_g(1)) > tol
    Ipb(:,:,2) = 1 - Ipb(:,:,2);
end   
if min_b(1) < 0 && abs(min_b(1)) > tol
    Ipb(:,:,3) = 1 - Ipb(:,:,3);
end
param.min_r = min_r; 
param.min_g = min_g; 
param.min_b = min_b; 
param.max_r = max_r; 
param.max_g = max_g; 
param.max_b = max_b; 
