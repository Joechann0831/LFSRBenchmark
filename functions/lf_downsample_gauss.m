function lf_ds = lf_downsample_gauss(lf_data_ori,scale)
% this function is used to downsample the lf data by scale factor of scale

[U,V,X,Y,C] = size(lf_data_ori);

lf_ds = zeros([U,V,(X-mod(X,scale))/scale,(Y-mod(Y,scale))/scale,C],'uint8');
sub_img_ds = zeros([(X-mod(X,scale))/scale,(Y-mod(Y,scale))/scale,C],'uint8');

for u = 1:U
    for v = 1:V
        sub_img = squeeze(lf_data_ori(u,v,:,:,:));
        sub_img = modcrop(sub_img,scale);
        sub_img_ds(:,:,1) = downsample_gauss(squeeze(sub_img(:,:,1)),scale);
        sub_img_ds(:,:,2) = downsample_gauss(squeeze(sub_img(:,:,2)),scale);
        sub_img_ds(:,:,3) = downsample_gauss(squeeze(sub_img(:,:,3)),scale);
        lf_ds(u,v,:,:,:) = sub_img_ds;
    end
end


end