function folder_list = get_folder_list(foldername)
% This function gives a list of folders withing the foldername.
% 
% Input: foldername: folder where the subfolders are stored
%
% Output: folder_list: cell_array containing a list of folders contained in 
%         foldername
% 
% Reuben Farrugia
% 11/8/2016
%

% Get information about the current foldername
info = dir(foldername);

% Get a cell list of files and folders contained within foldername
folder_list = {info(3:end).name}';
info(1:2) = [];
 
% Prune all the non-folders contained within info
folder_list([info(:).isdir]' == 0) = [];