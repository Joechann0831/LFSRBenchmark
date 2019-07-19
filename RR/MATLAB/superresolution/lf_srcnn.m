function SR_LF = lf_srcnn(LR_LF,mf)
addpath('MATLAB\matconvnet\');
addpath('MATLAB\matconvnet\matlab\');

%run matconvnet/matlab/vl_setupnn;
run vl_setupnn

model_filename = sprintf('DATA\\superresolution\\lfsrcnn\\lfsrcnn-model-x%d.mat',mf);

if ~exist(model_filename,'file')
    if mf == 2
        url =  'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=1boaQRfbQ7PYSrisHjgGxbdr_hDC_yN2b&export=download';
    elseif mf == 3
        url =  'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=10oM8wE-apPRzB7AyP6GRXFbDkn6aO0Sz&export=download';
    elseif mf == 4
        url =  'https://drive.google.com/a/um.edu.mt/uc?authuser=0&id=1U6QA_aRIo_wKr8Xu5Bzv7wxrIbF3At9l&export=download';
    end    
    % Download the required pca basis
    urlwrite(url,model_filename);
end
% Load the model filename
load(model_filename);

% Test that the convolution function works properly
try
    vl_nnconv(single(1),single(1),[]) ;
catch
    warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
    vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
        'enableImreadJpeg', false) ;
end

% Initialize the super-resolved light field
SR_LF = zeros(size(LR_LF),'uint8');

for u = 1:size(LR_LF,4)
    for v = 1:size(LR_LF,5)
        % Extract the current sub-aperture image
        Ilr = LR_LF(:,:,:,u,v);
        % Super-resolve the current sub-aperture image
        Isr = srcnn(Ilr,net);
        
        % Put the restored sub-aperture image in the restored light field
        SR_LF(:,:,:,u,v) = Isr;
    end
end

