function lf_gray = lf_rgb2gray(lf_data)
%this function is used to change the rgb ligth field to gray scale light
%field images

[U,V,X,Y,~] = size(lf_data);
lf_gray = zeros([U,V,X,Y],'uint8');

for u = 1:U
    for v = 1:V
        rgb_img = squeeze(lf_data(u,v,:,:,:));
        rgb_ycbcr = rgb2ycbcr(rgb_img);
        gray_img = squeeze(rgb_ycbcr(:,:,1));
        lf_gray(u,v,:,:) = gray_img;
    end
end
            

end