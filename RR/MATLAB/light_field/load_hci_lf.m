function LF = load_hci_lf(dataset_foldername,lf_name)

% Derive the foldername of the lf_name
lf_foldername = [dataset_foldername,lf_name];

for i = 0:80
    % Derive the image filename
    img_filename = sprintf('%s/input_Cam%.3d.png',lf_foldername,i);
    
    % Load the image
    I = imread(img_filename);
    
    % Determine the u and v parameters
    v = floor(i/9) + 1;
    u = mod(i,9)   + 1;
        
    % Store the light field LF(v,u,y,x,c)
    LF(v,u,:,:,:) = I;
end


