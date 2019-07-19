function Y = LF_bm3d_denoising(X)
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
addpath('MATLAB/noise_level_estimation/');

for i = 1:size(X,4)
    for j = 1:size(X,5)
        % Extract the noisy image z
        z = double(X(:,:,:,i,j))/255;
        
        % Estimate the noise variance
        nsig = NoiseLevel(255*z);
        
        nsig = mean(nsig);

        % Derive the estimate y
        [~,y] = CBM3D(1,z,nsig);
        
        % Put it in the right dynamic range
        y = uint8(y*255);
        
        % Put the restored sub-aperture image
        Y(:,:,:,i,j) = y;
    end
end
Y = uint8(Y);
