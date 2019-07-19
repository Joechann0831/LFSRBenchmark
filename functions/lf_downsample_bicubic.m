function lf_ds = lf_downsample_bicubic(lf_data_ori,scale)
% this function is used to downsample the lf data spatially

[U,V,X,Y,C] = size(lf_data_ori);

lf_ds = zeros([U,V,(X-mod(X,scale))/scale,(Y-mod(Y,scale))/scale,C],'uint8');

for u = 1:U
    for v = 1:V
        sub_img = squeeze(lf_data_ori(u,v,:,:,:));
        sub_img = modcrop(sub_img,scale);
        sub_img_ds = imresize(sub_img,1/scale);
        % downsample the lf data
        lf_ds(u,v,:,:,:) = sub_img_ds;
    end
end


end