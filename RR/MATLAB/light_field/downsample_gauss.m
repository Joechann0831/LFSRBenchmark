function img_ds = downsample_gauss(img_ori,scale)

%gaussion smoothing filtering

%img_smooth = imgaussfilt(img_ori,2);
H = fspecial('gaussian',[scale,scale],2);
img_smooth = imfilter(img_ori,H);

[X,Y] = size(img_smooth);
if scale == 3
    rad = floor(scale/2);
    img_ds = img_smooth((scale-rad):scale:(X-rad),(scale-rad):scale:(Y-rad));
else
    img_ds = zeros([X/2,Y/2],'uint8');
    for x = 1:X/2
        for y = 1:Y/2
            img_ds(x,y) = mean(mean(img_smooth((x-1)*2+1:x*2,(y-1)*2+1:y*2)));
        end
    end
end
% img_ds = img_smooth(3:5:X-2,3:5:Y-2);
end