function [B_r, B_g, B_b] = principal_basis_decoding(Ipb, param)

tol = 1E-5;

min_r = param.min_r; 
min_g = param.min_g; 
min_b = param.min_b; 
max_r = param.max_r; 
max_g = param.max_g; 
max_b = param.max_b; 

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
% Decompose the principal basis in rgb
B_r = Ipb(:,:,1);
B_g = Ipb(:,:,2);
B_b = Ipb(:,:,3);

% Reshape the RGB components of the principal basis as column vectors
B_r = B_r(:);
B_g = B_g(:);
B_b = B_b(:);

% Renormalize the B component
B_r = B_r.*(max_r - min_r) + min_r;
B_g = B_g.*(max_g - min_g) + min_g;
B_b = B_b.*(max_b - min_b) + min_b;



