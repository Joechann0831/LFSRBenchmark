%ok, GC for depth regularization

%class = 1xN vector of initial estimates

%unary = depthlabels x N of potential terms.(robust norm diff from initial)

%pairwise = sparse NxN link cost, positive for links, confidence & edge

%labelcost = depthlabels x depthlabels , cost of adj. depth diff. robust

%expansion = 1

 

%code assumes min_depth is 1

 

% load data

 

% image = pinhole

% depth = depth (1-255)

% confidence 

% all double. (note, im2double for depth does not work. use double()

pinhole_path        = '0-pinhole.png';
confidence_path     = '2-shear_depth_confidence.png';
depth_path          = '1-shear_depth_estimate.png';
confidence_path     = '4-corre_depth_confidence.png';
depth_path          = '3-corre_depth_estimate.png';

image = double(imread(pinhole_path));
depth = double(imread(depth_path));
confidence = double(imread(confidence_path)).^0.1;

 

min_depth = 1;

max_depth = 255;

 

[height width b] = size(depth);

 

class = depth(:)';

numnodes = size(class,2);

 

numdepth= max_depth - min_depth + 1;

 

unary = zeros(numdepth,numnodes);

 

flat_clip = 25;

for i = min_depth:max_depth

   

    unary(i,:) = min(abs(i - double(class)),flat_clip);

    

end

 

labelcost = zeros(numdepth,numdepth);

 

for i = min_depth:max_depth

    for j = min_depth:max_depth

   

        labelcost(i,j) = min(abs(i-j),flat_clip);

        

    end

end

 

%pairwise. goes down column first, width-1 * height + height-1 * width

 

 h = fspecial('sobel');

 

 image = im2double(image);

 

 output(:,:,1) = conv2(image(:,:,1),h,'same');

 output(:,:,2) = conv2(image(:,:,2),h,'same');

 output(:,:,3) = conv2(image(:,:,3),h,'same');

 

 vertical = abs(output(:,:,1)) + abs(output(:,:,2)) + abs(output(:,:,3));

 

 output(:,:,1) = conv2(image(:,:,1),h','same');

 output(:,:,2) = conv2(image(:,:,2),h','same');

 output(:,:,3) = conv2(image(:,:,3),h','same');

 

 horizontal = abs(output(:,:,1)) + abs(output(:,:,2)) + abs(output(:,:,3));

 

 %the confidence is quite crappy now. maybe i shouldn't use it.

 

 numentries = (width-1)*height + (height-1)*width;

 numentries = numentries*2 ;% make it sym

 

 %let ii = 

 ii = zeros(numentries,1);

 jj = zeros(numentries,1);

 ss = zeros(numentries,1);

 

 con_bias = 1;

 con_bias = 1;

 

 count = 1;

 for col = 1:width

     for row = 1:height-1

         

         total_edgestr = vertical(row,col) + vertical(row+1,col) + 0.2;

         local_con = confidence(row,col) + confidence(row+1,col);

         %the more confident it is, the more weight on the edge

         %if it is a strong edge in image, lower it

         local_weight = (local_con + con_bias)/total_edgestr;

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

         

         total_edgestr = vertical(row,col) + vertical(row,col+1) + 0.2;

         local_con = confidence(row,col) + confidence(row,col+1);

         %the more confident it is, the more weight on the edge

         %if it is a strong edge in image, lower it

         local_weight = (local_con + con_bias)/total_edgestr;

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

 

 pairwise = sparse(ii,jj,ss);

 

 

 [labels E Eafter] = GCMex(class, single(unary), pairwise, single(labelcost),1);

 

 results = reshape(labels,height,width);