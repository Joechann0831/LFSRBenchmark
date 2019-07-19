function imgH2 = iteration_once(imgH1,imgL0,scale,ds_flag)
% this function is used for once iteration 
% the downsampling method is an important factor that matters a lot
% so if the ds_flag = 1, we choose bicubic downsample
% otherwise we choose gauss downsample.

if ds_flag == 1
    imgL1 = imresize(imgH1,1/scale);
else
    imgL1 = downsample_gauss(imgH1,scale);
end

d1 = imgL0 - imgL1;
D1 = imresize(d1,scale);

imgH2 = imgH1+D1;

end