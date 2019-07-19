function SR_LF = lf_vdsr(LR_LF,mf)
addpath('MATLAB\matconvnet\');
addpath('MATLAB\matconvnet\matlab\');

%run matconvnet/matlab/vl_setupnn;
run vl_setupnn

addpath('MATLAB\superresolution\vdsr\utils\');
load('MATLAB\superresolution\vdsr\VDSR_Official.mat');

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
        Isr = vdsr(Ilr,model,mf);
        
        % Put the restored sub-aperture image in the restored light field
        SR_LF(:,:,:,u,v) = Isr;
    end
end

