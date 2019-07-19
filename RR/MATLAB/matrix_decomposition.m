function [B,C] = matrix_decomposition(X,k)

% Convert the 4D light field to a 3D tensor
X = reshape(X,[size(X,1),size(X,2),size(X,3)*size(X,4)]);

if nargin == 1
    % Determine the number of sub-aperture images
    k = size(X,3);
end
% Convert the 3D tensor to a matrix
X = reshape(X,[size(X,1)*size(X,2),size(X,3)]);
% Convert type of X to double
X = double(X);
% Compute singular value decomposition
[U,S,V] = svds(X,k);
% Derive the basis matrix B
B = U*S;
% Derive the coefficient matrix
C = V';

