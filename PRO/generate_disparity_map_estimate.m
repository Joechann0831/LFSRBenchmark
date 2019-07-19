function disparity_map = generate_disparity_map_estimate(disparity , U , V)
% This function is used to generate disparity map over the whole light
% field using the disparity map estimated for central view image.

[X,Y] = size(disparity);

cv_u = ceil(U/2);
cv_v = ceil(V/2);
disparity_map = zeros([U,V,X,Y,2]);

for u = 1:U
    for v = 1:V
        disparity_x = (u-cv_u)*disparity;% relative coordinates
        disparity_y = (v-cv_v)*disparity;
        disparity_map(u,v,:,:,1) = disparity_x;
        disparity_map(u,v,:,:,2) = disparity_y;
    end
end

        


end