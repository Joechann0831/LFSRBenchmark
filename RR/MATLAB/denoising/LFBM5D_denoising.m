function Y = LFBM5D_denoising(X,sig)

Y = zeros(size(X));

% Make sure that the library of bm3d is available
lfbm5d_foldername = 'MATLAB/denoising/LFBM5D/';

if ~exist(lfbm5d_foldername,'dir')
    if ~exist('MATLAB/denoising/','dir')
        mkdir('MATLAB/denoising/');
    end
    
    cd 'MATLAB/denoising/';
    
    % Clone the repository
    system('bash -c "git clone https://github.com/malain35/LFBM5D"');
    % This was tested on windows 10 using the linux subsystem.
    
    cd 'LFBM5d/';
    
    % Build the provided source
    system('bash -c "make [OMP=1]"');
    
    % Go back to where it was
    cd ..
    cd ..
    cd ..
end

% Create a noisy foldername
noisy_foldername = 'MATLAB/denoising/LFBM5D/LFNoisyDir/';
if ~exist(noisy_foldername,'dir')
    % Make sure the directory exists
    mkdir(noisy_foldername);
else
    delete_all_files(noisy_foldername);
end

% Put the distorted light field in this directory
lfbm5d_write_LF(noisy_foldername,X);

U = size(X,5); V = size(X,4);

if ~exist('MATLAB/denoising/LFBM5D/basicLF/','dir')
    mkdir('MATLAB/denoising/LFBM5D/basicLF/');
end

if ~exist('MATLAB/denoising/LFBM5D/denoisedLF/','dir')
    mkdir('MATLAB/denoising/LFBM5D/denoisedLF/');
end

% Derive the string command to be provided to call bash
command = sprintf('bash -c "MATLAB/denoising/LFBM5D/LFBM5Ddenoising none %d %d 1 1 row %d 2.7 MATLAB/denoising/LFBM5D/LFNoisyDir/ MATLAB/denoising/LFBM5D/basicLF/ MATLAB/denoising/LFBM5D/denoisedLF/ none 8 18 6 16 4 id sadct haar 0 16 18 6 8 4 dct sadct haar 0 opp 0 MATLAB/denoising/LFBM5D/outputMeasuresLFBM5D.txt"\n',U,V,sig);

% Launch the command line to compute the denoising
system(command);

% Read the light field which was restored using LFBM5D
Y = lfbm5d_read_LF('MATLAB/denoising/LFBM5D/denoisedLF/',U,V);

delete_all_files('MATLAB/denoising/LFBM5D/basicLF/');
delete_all_files('MATLAB/denoising/LFBM5D/denoisedLF/');

function delete_all_files(noisy_foldername)
% Clear all files from the file if it exists
dinfo = dir(noisy_foldername);
dinfo([dinfo.isdir]) = [];   %skip directories
if ~isempty(dinfo)
    filenames = fullfile(noisy_foldername, {dinfo.name});
    delete( filenames{:} );
end

function Y = lfbm5d_read_LF(foldername,U,V)

for v = 1:V
    for u = 1:U
        % Derive the filename
        filename = sprintf('%sSAI_%.2d_%.2d.png',foldername,v-1,u-1);
        I = imread(filename);
        Y(:,:,:,v,u) = I;
    end
end
Y  = uint8(Y);

function lfbm5d_write_LF(foldername,X)
% This writes the noisy light field
for v = 1:size(X,4)
    for u = 1:size(X,5)
        % Get the sub-aperture image
        I = X(:,:,:,v,u);
        % Derive the output filename
        filename = sprintf('%sSAI_%.2d_%.2d.png',foldername,v-1,u-1);
        % Write the sub-aperture images in the folder
        imwrite(I,filename);
    end
end


