function I_SR = srcnn(I_LR,net)
% This function is a wrapper script that allows us to use the code provided
% by the authors of the paper Chao Dong, Chen Change Loy, Kaiming He, Xiaoou 
% Tang. Image Super-Resolution Using Deep Convolutional Networks, IEEE 
% Transactions on Pattern Analysis and Machine Intelligence (TPAMI), 2015.
% The SRCNN was also used in a recent paper in lightfield super-resolution
% as mentioned by the same authors in Y. Yoon, H. G. Jeon, D. Yoo, J. Y. Lee
% and I. S. Kweon, "Learning a Deep Convolutional Network for Light-Field 
% Image Super-Resolution," 2015 IEEE International Conference on Computer 
% Vision Workshop (ICCVW), Santiago, 2015, pp. 57-65.
% This method receives the low-resolution lightfield LR_LF, and computes
% super-resolution so that each sub-aparture images has a dimension sz with
% a magnification factor provided in scale.
%
% Input: LR_LF - low-resolution lightfield
%        sz    - dimensions of each sub-aparture image
%        scale - magnification factor
%
% Output: SR_LF - the super-resolved lightfield image.
%
% Reuben Farrugia
% 30/07/2016

% Convert the low-resolution image in YCbCr color model
I_LR = rgb2ycbcr(I_LR);

% Initialize the super-resolved image
I_SR = I_LR;

% Extract the luminance component    
I_LR = I_LR(:,:,1);

% Convert the low-resolution image to double
I_LR = double(I_LR)/255;
    
% Compute the super-resolution
I_SR(:,:,1) = uint8(round(255*SRCNN_Matconvnet(I_LR,net)));

% Convert image in RGB
I_SR = ycbcr2rgb(I_SR);

function im_h_y = SRCNN_Matconvnet(im_l_y,net)

im_y = single(im_l_y);
convfea = vl_nnconv(im_y,net.layers{1}.weights{1},net.layers{1}.weights{2},'Pad',4);
convfea = vl_nnrelu(convfea);
convfea = vl_nnconv(convfea,net.layers{3}.weights{1},net.layers{3}.weights{2},'Pad',2);
convfea = vl_nnrelu(convfea);
im_h_y  = vl_nnconv(convfea,net.layers{5}.weights{1},net.layers{5}.weights{2},'Pad',2);

