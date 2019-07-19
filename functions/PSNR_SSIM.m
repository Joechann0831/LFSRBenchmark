function [average_PSNR,PSNR_var,average_SSIM,SSIM_var,central_PSNR,central_SSIM,PSNRs] = PSNR_SSIM(lf_gt,lf_estimate,scale)

[U,V,~,~] = size(lf_estimate);
PSNRs = zeros([U,V]);
SSIMs = zeros([U,V]);

for u = 1:U
    for v = 1:V
        sub_gt = squeeze(lf_gt(u,v,:,:));
        sub_gt = modcrop(sub_gt,scale);
        PSNRs(u,v) = psnr(sub_gt,squeeze(lf_estimate(u,v,:,:)));
        SSIMs(u,v) = ssim(sub_gt,squeeze(lf_estimate(u,v,:,:)));
    end
end

average_PSNR = mean(PSNRs(:));
average_SSIM = mean(SSIMs(:));
PSNR_var = var(PSNRs(:));
SSIM_var = var(SSIMs(:));
central_PSNR = PSNRs(5,5);
central_SSIM = SSIMs(5,5);


end