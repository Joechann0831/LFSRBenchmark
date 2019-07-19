function [D_lrlf, D_hrlf] = load_coupled_dictionaries( param, out_filename)
% This script loads the coupled dictionaries D_lrlf and D_hrlf. In summary
% we first load the mat file filename which contains a dict structure
% containing 2000 patch volumes per lightfield. We then iterature through
% all 9 lightfields considered in this experiment for the HCI lightfields
% and only the patches which do not correspond to patchs in lf_name are
% included in the coupled dictionaries
%
% Input: filename - filename of the mat file containing the patches
%        lf_name  - name of the lightfield being super-resolved
%
% Output: D_lrlf - low-resolution lightfield dictionary
%         D_hrlf - high-resolution lightfield dictionary
%
% Reuben Farrugia
% Date: 29/07/2016

if ~exist(out_filename,'file')
    %fprintf('Loading the coupled dictionary from %s.\n', filename);
    list_filename = 'DATA/train_test_split.mat';

    list = load(list_filename);

    list = list.train_list;
    Npts = 1000;

    % Generate a mat file containing te required dictionaries
    [D_hrlf,D_lrlf] = generate_coupled_dictionary(param,list,Npts);
    
    %--------------------------------------------------------------------------
    % This part constructs the dictionary. The dict contains a number of patch
    % volues (2000 for each lightfield). However, in this simulation we
    % construct coupled dictionaries using the leave one out methodologies
    % where only patches that do not correspond to patches from the lf_name
    % lightfield are included in the dictionary. 
%     D_hrlf = []; D_lrlf = [];
%     for n = 1:N
%         Dh = dict{n}.Dlf_H; Dl = dict{n}.Dlf_L;
% 
%         % Concatenate the dictionaries
%         D_hrlf = [D_hrlf, Dh]; D_lrlf = [D_lrlf, Dl];
%     end
    fprintf('Save the dictionaries\n..');
    save(out_filename,'D_hrlf','D_lrlf','-v7.3');
else
    load(out_filename);
end

