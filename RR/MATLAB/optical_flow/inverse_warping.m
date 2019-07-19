function A = inverse_warping(X,u,v)
% This function will align all the sub-aperture images to the centre view
% using the flow vectors (u,v) computed using optical flow
%
% Reuben Farrugia
% 26/1/2018
%

% Initialize the aligned light field
A = zeros(size(X),'uint8');

k = 0;
for i = 1:size(X,4)
    for j = 1:size(X,5)
        %--- Display the progress of the sift-flow computation
        msg = sprintf('  Light field alignment: %6.2f%%', k/(size(X,4)*size(X,5))*100);
        fprintf('%s',msg);
        lengthLastMsg = length(msg);
        pause(0.005);
        
        % Extract the non-aligned view (i,j)
        It = X(:,:,:,i,j);
        
        % Warp pixels from the sournce It to theIt aligned light field A
        for x = 1:size(It,2)
            for y = 1:size(It,1)
                % Derive the new values of x and y
                y_new = min(max(round(y - v(y,x,i,j)),1),size(It,1));
                x_new = min(max(round(x - u(y,x,i,j)),1),size(It,2));
                
                % Compute the warping
                A(y,x,:,i,j) = It(y_new,x_new,:);
            end
        end
        
        % Increment counter
        k = k + 1;
        
        %--- Clear the last entry
        fprintf(repmat('\b', 1, lengthLastMsg));
    end
end