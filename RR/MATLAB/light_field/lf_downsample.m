function LR_LF = lf_downsample(HR_LF,mf,flag)
% This function is used to downsample the light field using a scale 1/mf
% and then rescale it using bicubic interpolation.

if nargin == 2
    flag = 0;
end

for u = 1:size(HR_LF,4)
    for v = 1:size(HR_LF,4)
        % Get the current view
        I = HR_LF(:,:,:,u,v);
        % Downsample I by a scale 1/mf
        I = imresize(I,1/mf);
        if flag == 0
            % Rescale I using bicubic interpolation
            Ilr = imresize(I,mf);
            % Crop the interpolated image
            Ilr = Ilr(1:size(HR_LF,1),1:size(HR_LF,2),:);
            % This was giving a pixelwise misalignment when the resolution
            % was not divisible by the magnification factor
            %LR_LF(:,:,:,u,v) = imresize(I,[size(HR_LF,1),size(HR_LF,2)]);
            % Put the degraded light field into the tensor
            LR_LF(:,:,:,u,v) = Ilr;
        else
            LR_LF(:,:,:,u,v) = I;
        end
    end
end