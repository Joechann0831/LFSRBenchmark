function    [V,D,M] = pca(X,N)
% This script computes the pca on the input dictionary X. In this work we
% developed our own implementation of PCA to avoid licensing issues since
% most of the time no statistical toolbox license was available. It returns
% the eigenvectors V, the eigenvalues D and the mean vector M.

% Make sure the input vector is single precision
X = single(X);

% Compute the mean vector
M = mean(X,2);

% Center the data
X = X - repmat(M,[1,size(X,2)]);

% Compute the covariance matrix
C = double(X*X')';
    
% Derive the eigen values and eigenvectors of the covariance matrix
[V,D] = eigs(C,N); 