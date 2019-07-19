function lf_ds = lf_downsample_bicubic(lf_data_ori,scale)
% this function is used to downsample the lf data spatially and upsample it
% back to the original resolution after modcrop

[U,V,X,Y,C] = size(lf_data_ori);

lf_ds = zeros([U,V,(X-mod(X,scale)),(Y-mod(Y,scale)),C],'uint8');

for u = 1:U
    for v = 1:V
%         if scale == 3
%             sub_img_ds = uint8(img_downsample_3(squeeze(lf_data_ori(u,v,:,:,:))));
%         else
%             sub_img_ds = imresize(squeeze(lf_data_ori(u,v,:,:,:)),1/scale);
%         end
        sub_img = squeeze(lf_data_ori(u,v,:,:,:));
        sub_img = modcrop(sub_img,scale);
        sub_img_ds = imresize(sub_img,1/scale);
        %downsample the lf data
        lf_ds(u,v,:,:,:) = imresize(sub_img_ds,scale);
    end
end


end