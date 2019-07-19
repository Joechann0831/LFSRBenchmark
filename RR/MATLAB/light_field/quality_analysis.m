function [psnr_v,ssim_v] = quality_analysis(HR_LF, SR_LF)
% This script is a simple script which evaluates the quality of the
% super-resolved lightfield. This function returns two quality metrics PSNR
% and SSIM

%fprintf('Computing the quality metrics\n');

% Determine the number of sub-aparture images
N = size(HR_LF,4)*size(HR_LF,5);

% Initialize the quality arrays to contain a quality per sub-aparutre image
psnr_v = zeros(N,1); 
ssim_v = zeros(N,1);

HR_LF = uint8(HR_LF); SR_LF = uint8(SR_LF);
k = 1;
for n1 = 1:size(HR_LF,4)
    for n2 = 1:size(HR_LF,5)
        if n1 >= 3 && n1 <= 7 && n2 >= 3 && n2 <=7
            % Load the high-resolution sub-aparture image
            I1 = rgb2gray(HR_LF(:,:,:,n1,n2));
            % Load the super-resolved sub-aparture image
            I2 = rgb2gray(SR_LF(:,:,:,n1,n2));
            
            % Crop the images to ignore the borde
            I1 = I1(21:end-20, 21:end-20);
            I2 = I2(21:end-20, 21:end-20);
            
            % Compute the ssim quality metric
            psnr_v(k) = psnr(I1,I2);
            ssim_v(k) = ssim(I1,I2);
        end
        k = k + 1;
    end
end

psnr_v(psnr_v == 0) = [];
ssim_v(ssim_v == 0) = [];
% Average the quality metrics to get a unique metric for all sub-aparture
% images
psnr_v = mean(psnr_v);
ssim_v = mean(ssim_v);