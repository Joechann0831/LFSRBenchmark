function [lf_names, datasets] = read_configuration(cfg_filename)
% This file is used to read information contained in the configuration
% file. The data is stored in the format lf_names, datasets

% Open the file to be read from 
fid = fopen(cfg_filename,'r');

k = 1;
while ~feof(fid)
    % Read a line
    line = fgetl(fid);
    % Derive the index of commas
    dlm_idx = strfind(line,',');
    % Extract the light field name
    lf_names{k} = strtrim(line(1:dlm_idx-1));
    % Extract the dataset name
    datasets{k} = strtrim(line(dlm_idx+1:end));
    k = k + 1;
end


