function [X_y,X_cb,X_cr] = LF_rgb2ycbcr(X)
% The LF_rgb2ycbcr script is used to convert each sub-aperture image in X
% into the YCbCr color channel. Each color channel is given separately as
% X_y for luminance and X_cb,X_cr for the chrominance component
%
% Reuben Farrugia
% 26/1/2018
%
%----- Convert image in YCbCr color model
X_y = zeros(size(X,1),size(X,2),size(X,4)*size(X,5),'uint8');
X_cb = X_y;
X_cr = X_y;

% Convert the light field in the YCbCr color model
k = 1;
for i = 1:size(X,4)
    for j = 1:size(X,5)
        % Get the sub-aperture image in rgb color model
        I_rgb = uint8(reshape(X(:,:,:,i,j),[size(X,1),size(X,2),size(X,3)]));
        % Derive the I_ycbcr image
        I_ycbcr = rgb2ycbcr(I_rgb);
        % Convert the RGB image into Ycbcr
        % Put the luma component in the light-field
        X_y(:,:,k) = I_ycbcr(:,:,1); % Luminance component
        X_cb(:,:,k) = I_ycbcr(:,:,2);% First chroma component
        X_cr(:,:,k) = I_ycbcr(:,:,3);% Second chroma component
        % Increment the counter
        k = k + 1;
    end
end