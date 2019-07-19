function depth_output  = DEPTH_MRF2_s(depth_buffer, confi_buffer, occu,...
    IM_Pinhole, height, width, depth_resolution)

% class = 1xN vector of initial estimates
% unary = depthlabels x N of potential terms.(robust norm diff from initial)
% pairwise = sparse NxN link cost, positive for links, confidence & edge
% labelcost = depthlabels x depthlabels , cost of adj. depth diff. robust
% expansion = 1
tic;

img_size = width * height;
confi_buffer(~isfinite(confi_buffer)) = 0;
confi_buffer(depth_buffer < 0.1*depth_resolution | depth_buffer > 0.9*depth_resolution) = 0;

% confidence gaussian for UNARY
confi_gaus_sigma_mult = 0.1 * depth_resolution;
% labelcost gaussian for LABELCOST
label_gaus_radius = depth_resolution;
label_gaus_sigma = 0.05 * depth_resolution;

%% Class: Initial Depth Labels
CLASS = reshape(depth_buffer, 1, img_size);
%% Unary: Confidence
% confidence  
confi_val = 1 - confi_buffer;
confi_val = min(1, max(0.01, confi_val));

d = 1:depth_resolution;
buffer = depth_buffer(1:img_size)';
[D, BUFFER] = meshgrid(d, buffer);
confi_gaus_sigma = repmat(confi_gaus_sigma_mult*confi_val(:), [1 depth_resolution]);
UNARY = exp(-(D - BUFFER).^2./(2*confi_gaus_sigma.^2))'./(sqrt(2*pi)*confi_gaus_sigma');
UNARY = UNARY ./ repmat(max(UNARY), [depth_resolution 1]);
UNARY = 1 - UNARY;

%% PAIRWISE: Image pixel color smoothness enforcement
h = fspecial('sobel');

image = im2double(IM_Pinhole);
output(:,:,1) = conv2(image(:,:,1),h,'same');
output(:,:,2) = conv2(image(:,:,2),h,'same');
output(:,:,3) = conv2(image(:,:,3),h,'same');

vertical = abs(output(:,:,1)) + abs(output(:,:,2)) + abs(output(:,:,3));
output(:,:,1) = conv2(image(:,:,1),h','same');
output(:,:,2) = conv2(image(:,:,2),h','same');
output(:,:,3) = conv2(image(:,:,3),h','same');

horizontal = abs(output(:,:,1)) + abs(output(:,:,2)) + abs(output(:,:,3));
numentries = (width-1)*height + (height-1)*width;
numentries = numentries*2 ;% make it sym

ii = zeros(numentries,1);
jj = zeros(numentries,1);
ss = zeros(numentries,1);

edge_str = 0.2;
edge_padding = 0.1;

im_edge = edge(rgb2gray(IM_Pinhole), 'canny', 0.1, 1);
canny_weight = max(max(vertical+horizontal));
occu_weight = 1e5*max(max(vertical+horizontal));

count = 1;
for col = 1:width
     for row = 1:height-1
         total_edgestr = (abs(vertical(row, col) - vertical(row+1, col)) + ...
             canny_weight * abs(im_edge(row, col) - im_edge(row+1, col)) + ...
             occu_weight * abs(occu(row, col) - occu(row+1, col)) + edge_padding) * edge_str;         
         %the more confident it is, the more weight on the edge
         %if it is a strong edge in image, lower it
         local_weight = 1/total_edgestr;
         ss(count) = local_weight;
         ii(count) = row + (col-1)*height;
         jj(count) = row + 1 + (col-1)*height;
         count = count + 1;
         ss(count) = local_weight;
         jj(count) = row + (col-1)*height;
         ii(count) = row + 1 + (col-1)*height;
         count = count + 1;

     end
end

for col = 1:width-1
     for row = 1:height
         total_edgestr = (abs(horizontal(row, col) - horizontal(row, col+1)) + ...
            canny_weight * abs(im_edge(row, col) - im_edge(row, col+1)) + ...
            occu_weight * abs(occu(row, col) - occu(row, col+1)) + edge_padding) * edge_str;         
         local_weight = 1/total_edgestr;
         ss(count) = local_weight;
         ii(count) = row + (col-1)*height;
         jj(count) = row + (col)*height;
         count = count + 1;
         ss(count) = local_weight;
         jj(count) = row + (col-1)*height;
         ii(count) = row + (col)*height;
         count = count + 1;

     end
end

PAIRWISE = sparse(ii,jj,ss);

%% LABELCOST: Depth change smoothness enforcement
LABELCOST = eye(depth_resolution,depth_resolution)                        ;
h         = fspecial('gaussian',[label_gaus_radius 1],label_gaus_sigma)   ;
LABELCOST = imfilter(LABELCOST,h,'symmetric')                             ;
LABELCOST = LABELCOST ./ repmat(max(LABELCOST), [depth_resolution 1]);
LABELCOST = 1 - LABELCOST                                                 ;
LABELCOST = min(0.3, LABELCOST);

labels    = GCMex(CLASS, single(UNARY), PAIRWISE, single(LABELCOST),1)    ;

depth_output        = reshape(labels,height,width)                        ;

fprintf('Completed in %.3f seconds\n', toc);

