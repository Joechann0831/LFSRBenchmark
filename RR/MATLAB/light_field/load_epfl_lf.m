function LF = load_epfl_lf(dataset_foldername,lf_name)

% Derive the foldername of the lf_name
lf_filename = sprintf('%s%s.mat',dataset_foldername,lf_name);

% Load the matlab filename
data = load(lf_filename);

% Extract the light-field image
LF = data.LF;

LF = uint8(round(single(LF(3:11,3:11,:,:,1:3))/(2^16-1)*(2^8-1)));

% Clear the data structure
clearvars data;

