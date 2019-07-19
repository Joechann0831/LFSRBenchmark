function pca_basis = learn_pca_basis(D_lrlf, D_hrlf, param)
% This script is used to learn the pca basis for the coupled dictionaries.
% The pca basis are used to project the low-resolution samples on the
% low-resolution subspace and the high-resolution samples are projected on
% the high-resolution subspace. The dimension of the sub-space is
% significantly lower than the dimension of the patch volume. However, here
% we are returning a sub-space with dimension equivalent to the dimension
% of the patch volume i.e. 8100 for a patch-size of 10 and 81 sub-aparture
% images. The eigenvalues for each subpsace are returned into the pca_basis
% structure since the dimensionality reduction is computed not in this
% script. Therefore, if one needs to reduce the dimension of one of the
% subspaces you have to consider the magnitude of the eigen-values and
% choose the first N eigenvectors of the pca basis.
% 
% Input: D_lrlf - low-resolution patch-volume dictionary
%        D_hrlf - high-resolution patch-volume dictionary
%        mat_filename - mat filename where the data will be stored
%
% Output: pca_basis - a structure containing the eigenvectors, eigenvalues,
% and mean vectors for each subspace. In detail
%
% pca_basis.El : Eigenvectors that can be used to project low-quality
%                samples onto the low-quality subspace
% pca_basis.Ml : Mean vector of the low-quality patch-volumes
% pca_basis.Dl : Eigenvalues of the low-quality subspace
% pca_basis.Eh : Eigenvectors that can be used to project high-quality
%                samples onto the high-quality subspace
% pca_basis.Mh : Mean vector of the high-quality patch-volume
% pca_basis.Dh : Eigenvalues of the high-quality subspaces
%
% Reuben Farrugia
% Date: 29/7/2016
%
% NOTE: This function will be slow when computed the first time for lf_name
% at a specified magnification factor since the learned pca basis will be
% stored in a mat file and loaded from that point onwards.

fprintf('Learning the PCA basis.\n');

% Compute pca dimensionality reduction using all points
[El, Dl, Ml] = pca(D_lrlf, param.Nl);

% Compute pca dimensionality redunction using all points
[Eh, Dh, Mh] = pca(D_hrlf, param.Nh);
    
% Put the basis and mean into a structre
pca_basis.El = El; pca_basis.Ml = Ml; pca_basis.Dl = diag(Dl);
pca_basis.Eh = Eh; pca_basis.Mh = Mh; pca_basis.Dh = diag(Dh);  
