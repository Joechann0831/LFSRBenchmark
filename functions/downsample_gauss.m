function img_ds = downsample_gauss(img_ori,scale)
% This function is used for gaussian downsampling with scale = 2 or 3
% With scale = 2 or 3, the kernel size 3 is an appropriate size

%gaussion smoothing filtering
H = fspecial('gaussian',[3,3],2);
img_smooth = imfilter(img_ori,H,'replicate');

[X,Y] = size(img_smooth);
if scale == 3
    rad = floor(scale/2);
    img_ds = img_smooth((scale-rad):scale:(X-rad),(scale-rad):scale:(Y-rad));
else
    img_smooth_ul = double(img_smooth(1:2:end,1:2:end));
    img_smooth_ur = double(img_smooth(1:2:end,2:2:end));
    img_smooth_bl = double(img_smooth(2:2:end,1:2:end));
    img_smooth_br = double(img_smooth(2:2:end,2:2:end));
    img_ds = uint8((img_smooth_ul + img_smooth_ur + img_smooth_bl + img_smooth_br) / 4.0);
end
end