function Y = LF_ycbcr2rgb(X_y,X_cb,X_cr)
% The LF_ycbcr2rgb script is used to convert each sub-aperture image in 
% into the RGB color channel. Each color channel is given separately as
% X_y for luminance and X_cb,X_cr for the chrominance component
%
% Reuben Farrugia
% 26/1/2018
%
N1 = sqrt(size(X_y,3));
N2 = sqrt(size(X_y,3));

% Initialize the light field in RGB space
Y = zeros(size(X_y,1),size(X_y,2),3,N1,N2,'uint8');

k = 1;
for i = 1:N1
    for j = 1:N2
        % Get the sub-aperture image in rgb color model
        I_ycbcr(:,:,1) = uint8(X_y(:,:,k));
        I_ycbcr(:,:,2) = uint8(X_cb(:,:,k));
        I_ycbcr(:,:,3) = uint8(X_cr(:,:,k));
        
        % Get the rgb image
        I_rgb = ycbcr2rgb(I_ycbcr);
        
        % Put the sub-aperture image in the light field
        Y(:,:,:,i,j) = I_rgb;

        % Increment counter
        k = k + 1;
    end
end

                               
