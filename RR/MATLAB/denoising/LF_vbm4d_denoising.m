function Y = LF_vbm4d_denoising(X)
% This function computed BM3D to restore the noise from the input light
% field
%

% Make sure that the library of bm3d is available
vbm4d_foldername = 'MATLAB/denoising/VBM4D/';

if ~exist(vbm4d_foldername,'dir')
    if ~exist('MATLAB/denoising/','dir')
        mkdir('MATLAB/denoising/');
    end
    % Set the url where the library is stored
    url = 'http://www.cs.tut.fi/~foi/GCF-BM3D/VBM4D_v1.zip';
    
    % Download the toolbox
    urlwrite(url,'MATLAB/denoising/VBM4D.zip');
    
    % Unip the file
    unzip('MATLAB/denoising/VBM4D.zip',vbm4d_foldername);
    
    % Delete the zipfile
    delete('MATLAB/denoising/VBM4D.zip');
end

% Add the bm3d library to the path
addpath('MATLAB/noise_level_estimation/');
addpath('MATLAB/denoising/VBM4D/');

profile = 'lc';      % V-BM4D parameter profile
                     %  'lc' --> low complexity
                     %  'np' --> normal profile
do_wiener = 1;       % Wiener filtering
                     %   1 --> enable Wiener filtering
                     %   0 --> disable Wiener filtering
sharpen = 1;         % Sharpening
                     %   1 --> disable sharpening
                     %  >1 --> enable sharpening
deflicker = 1;       % Deflickering
                     %   1 --> disable deflickering
                     %  <1 --> enable deflickering
verbose = 0;         % Verbose mode

% Get centre view
Ic = X(:,:,:,5,5);

% Estimate the noise variance
nsig = NoiseLevel(double(Ic));
        
nsig = mean(nsig);

% Reshape the lighfield to 4D
z = reshape(X,[size(X,1),size(X,2),size(X,3),size(X,4)*size(X,5)]);

% Scale the video
z = double(z)/255;

% Compute vbm4d
y = vbm4d( z , nsig/255, profile, do_wiener, sharpen, deflicker, verbose );

% Restore the 5D lightfield
Y = reshape(y,size(X));

% Convert to uint8
Y = uint8(Y*255);
