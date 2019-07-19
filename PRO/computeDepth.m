function depth_output = computeDepth(data , d_max)

% CONTACT:
% Ting-Chun Wang (tcwang0509@berkeley.edu)

% TERMS OF USE : 
% Any scientific work that makes use of our code should appropriately
% mention this in the text and cite our ICCV 2015 paper. For commercial
% use, please contact us.

% PAPER TO CITE:
% Ting-Chun Wang, Alexei A. Efros, and Ravi Ramamoorthi.
% Occlusion-aware depth estimation using light-field cameras. 
% In Proceedings of International Conference on Computer Vision (ICCV), 2015.

% BIBTEX TO CITE:
% @inproceedings{wang2015occlusion,
%   title={Occlusion-aware depth estimation using light-field cameras.},
%   author={Wang, Ting-Chun and Efros, Alexei and Ramamoorthi, Ravi},
%   booktitle={Proceedings of the IEEE International Conference on Computer Vision (ICCV)},
%   year={2015}
% }

% Copyright (c) 2015
% All rights reserved.
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution      
%     * Proper citation to the paper above
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

addpath required
addpath mex
addpath(genpath('regularize/'))

%% Parameters
% d_max = 1;          % maximum disparity between neighboring cameras
d_res = 200 * d_max + 1;        % depth resolution
v_max = 100;        % maximum correspondence cue response
f_max = 0.01;       % maximum refocus cue response
dilate_amount = 4;  % dilate amount on edge detection
occ_thre = 0.1;     % threshold to be regard as occlusion

% parameters from input
UV_diameter = size(data, 4);                           % angular resolution
UV_radius   = floor(UV_diameter/2);                    % half angular resolution
h           = size(data, 1);                           % spatial image height
w           = size(data, 2);                           % spatial image width
LF_y_size   = h * UV_diameter;                         % total image height
LF_x_size   = w * UV_diameter;                         % total image width
LF_Remap    = reshape(permute(data, ...
               [4 1 5 2 3]), [LF_y_size LF_x_size 3]); % the remap image the very big one it contains 9*9 images
im          = data(:,:,:,UV_radius+1,UV_radius+1);     % the pinhole image

%% Edge detection and orientation computation
tic; disp('1. Edge detection')
im_edge = edge(im(:,:,1), 'canny') | edge(im(:,:,2), 'canny') | edge(im(:,:,3), 'canny');
orient = skeletonOrientation(im_edge, [5 5]);
orient(~im_edge) = -100; 
dir = imdilate(orient, strel('disk', dilate_amount));
fprintf('Edge detecion completed in %.2f sec\n', toc);

%% Initial depth estimation
tic; disp('2. Initial depth estimation')
[depth, d_var, d_cost, d_foc] = ...
    occlusionDetection_mex(w, h, UV_diameter, LF_Remap, -d_max, d_max, d_res, dir*pi/180);
depth = double(depth);
fprintf('Initial depth estimation completed in %.2f sec\n', toc);

%% Occlusion cue computation
tic; disp('3. Occlusion cue computation')
occlusion = occCompute(d_var, d_foc, im_edge, dir, v_max, f_max, occ_thre);
fprintf('Occlusion cue computation completed in %.2f sec\n', toc);

%% Final depth estimation
tic; disp('4. Final depth estimation')
confidence   = confCompute(d_cost, h, w);
depth_output = DEPTH_MRF2_s(depth, confidence, occlusion, im, h, w, d_res);
depth_output = medfilt2(depth_output, [3 3]);
fprintf('Final depth estimation completed in %.2f sec\n', toc);

%% Save the data
disparity = - d_max + 0.01 * depth_output;
depth_output = disparity;
% disparity_reg=(disparity+1)/2;
% file_out = strcat('results\', name, '\dis_SR.mat');
% save(file_out, 'disparity', '-mat');
% imshow(disparity_reg);
% img_out = strcat('results\', name, '\dis_SR.png');
% imwrite(disparity_reg, img_out);
% disp('Storage Complete!');