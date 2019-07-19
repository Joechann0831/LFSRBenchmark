clc; close all; clear all;
% This script computes batch processing of the application to be specified.
% application defines the application to be processed.
% application: SR - super-resolution

application = 'SR';

if strcmp(application,'SR')
    % Specify the magnification factors to consider
    mfs = [2,4];
    % Specify the algorithms to be compared
    methods = {'pb-vdsr'};
    
    for i = 1:size(mfs,2)
        for j = 1:size(methods,2)
            mf = mfs(i);
            method = methods{j};
            % Compute the light field super-resolution analysis
            LF_super_resolution_analysis(method,mf,true);
        end
    end
end