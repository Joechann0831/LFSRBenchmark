function create_LF_video()
clc; close all; clear all;

%--------------------------------------------------------------------------
% Configuration
%--------------------------------------------------------------------------
mf = 3;
method = 'pca_rr';
lf_name = 'Fountain_&_Bench';

% Derive the light field filename
LF_filename = sprintf('../RESULTS/superresolution/x%d/LF/%s/%s.mat',mf,method,lf_name);

% Load the data
data = load(LF_filename); SR_LF = data.SR_LF; clearvars data;

% Reshape SR_LF into a sequence of video frames
SR_LF = reshape(SR_LF,[size(SR_LF,1),size(SR_LF,2),size(SR_LF,3),size(SR_LF,4)*size(SR_LF,5)]);

% Derive the output video filename
out_video_filename = sprintf('../RESULTS/superresolution/x%d/LF/%s/%s.avi',mf,method,lf_name);
v = VideoWriter(out_video_filename);
v.FrameRate = 5;

implay(SR_LF);

open(v);
writeVideo(v,SR_LF)
close(v);