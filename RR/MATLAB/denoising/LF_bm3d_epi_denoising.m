function Y = LF_bm3d_epi_denoising(X,sig)
% This function computed BM3D to restore the noise from the input light
% field
%

Y = zeros(size(X));
% Make sure that the library of bm3d is available
bm3d_foldername = 'MATLAB/denoising/BM3D/';

if ~exist(bm3d_foldername,'dir')
    if ~exist('MATLAB/denoising/','dir')
        mkdir('MATLAB/denoising/');
    end
    % Set the url where the library is stored
    url = 'http://www.cs.tut.fi/~foi/GCF-BM3D/BM3D.zip';
    
    % Download the toolbox
    urlwrite(url,'MATLAB/denoising/BM3D.zip');
    
    % Unip the file
    unzip('MATLAB/denoising/BM3D.zip',bm3d_foldername);
    
    % Delete the zipfile
    delete('MATLAB/denoising/BM3D.zip');
end

% Add the bm3d library to the path
addpath('MATLAB/denoising/BM3D/');

% Determine the number of sai vertically
N = size(X,4);
for i = 1:N
    for y = 1:size(X,1)
        % Form a horizontal epi
        epi = permute(reshape(X(y,:,:,i,:),[size(X,2),size(X,3),size(X,5)]),[3,1,2]);
        % Clean the epi using bm3d
        z = double(epi)/255;
        [~,epi] = CBM3D(1,z,sig);
        % Reshape the epi
        epi = permute(epi,[2,3,1])*255;
        Y(y,:,:,i,:) = epi;
    end
end
N = size(X,5);
for i = 1:N
    for x = 1:size(X,2)
        % Form a vertical epi
        epi = permute(reshape(Y(:,x,:,:,i),[size(X,1),size(X,3),size(X,4)]),[3,1,2]);
        % Clean the epi using bm3d
        z = double(epi)/255;
        [~,epi] = CBM3D(1,z,sig);
        % Reshape the epi
        epi = permute(epi,[2,3,1])*255;
        Y(:,x,:,:,i) = epi;
    end
end
Y = uint8(Y);
