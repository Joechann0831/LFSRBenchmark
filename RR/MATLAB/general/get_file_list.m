function list = get_file_list(foldername,ext)
% The get_file_list derives a cell list of files contained into the folder
% foldername. Only the files having the extension indicated in ext string
% are included in the list of cells.
%
% Reuben Farrugia
% Data: 29/7/2016

% Derive the information about this folder 
info = dir([foldername,'*.',ext]);
% Return the list of images with extension ext
list = {info(:).name}';