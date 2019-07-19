function lf_gray_ds = lf_gray_bicubic_ds(lf_gray,scale)

[U,V,X,Y] = size(lf_gray);

lf_gray_ds = zeros([U,V,(X-mod(X,scale))/scale,(Y-mod(Y,scale))/scale],'uint8');

for u = 1:U
    for v = 1:V
        sub_img = squeeze(lf_gray(u,v,:,:));
        sub_img = modcrop(sub_img,scale);
        sub_img_ds = imresize(sub_img,1/scale);
        lf_gray_ds(u,v,:,:) = sub_img_ds;
    end
end

end