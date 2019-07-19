function LF = load_stanford_lf(dataset_foldername,lf_name)
% This function is used to lead the stanford light fields (v,u,y,x,c)

% Derive the foldername of the lf_name
lf_foldername = [dataset_foldername,lf_name];

for v = 1:9
    for u = 1:9
        % Derive the image filename
        img_filename = sprintf('%s/IMG_%d_%d.png',lf_foldername,v,u);
        
        % Load the image
        I = imread(img_filename);
        
        u_new = u;
        v_new = 9 - v + 1;
        
        % Store the light field LF(v,u,y,x,c)
        LF(v_new,u_new,:,:,:) = I;
    end
end