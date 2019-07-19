function [u_,v_] = forward_optical_flow(LR_LF)
% This function uses the SIFTflow optical flow library to compute the flow
% vector for each sub-aperture image with respect to the center view. It
% takes the LR_LF as an input i.e. the non-aligned light field and returns
% the flow vectors u_ and v_.
%
% Reuben Farrugia
% 26/1/2018
%

% Derive the foldername where siftflow will be stored
siftflow_foldername = 'MATLAB/optical_flow/SIFTflow/';

if ~exist(siftflow_foldername,'dir')
    
    % Derive the siftflow filename
    siftflow_filename = 'MATLAB/optical_flow/SIFTflow.zip';
    
    % Derive the url of siftflow
    url = 'https://people.csail.mit.edu/celiu/SIFTflow/SIFTflow.zip';
    
    % Download the siftflow code from web
    urlwrite(url,siftflow_filename);
    
    % Unzip the folder
    unzip(siftflow_filename,'MATLAB/optical_flow/');
    
    % Delete the zip file for cleaning
   delete(siftflow_filename);
end

% Add the path to SIFTflow and other functions
addpath('MATLAB/optical_flow/SIFTflow/');
addpath('MATLAB/optical_flow/SIFTflow/mexDenseSIFT/');
addpath('MATLAB/optical_flow/SIFTflow/mexDiscreteFlow/');

%--- Parameter configuration of SIFTflow
SIFTflowpara.alpha=2*255;
SIFTflowpara.d=40*255;
SIFTflowpara.gamma=0.005*255;
SIFTflowpara.nlevels=4;
SIFTflowpara.wsize=2;
SIFTflowpara.topwsize=10;
SIFTflowpara.nTopIterations = 60;
SIFTflowpara.nIterations= 30;
%---

% Initialize the flow vectors
u_ = zeros(size(LR_LF,1),size(LR_LF,2),size(LR_LF,4),size(LR_LF,5));
v_ = zeros(size(LR_LF,1),size(LR_LF,2),size(LR_LF,4),size(LR_LF,5));

% Parameter specification for SIFTflow (as per default)
cellsize=3;
gridspacing=1;

% Extract the center view
Ic = rgb2gray(LR_LF(:,:,:,5,5)); % Grayscale component only

% Compute the dense sift feature extraction
sift1 = mexDenseSIFT(Ic,cellsize,gridspacing);

k = 0;
for i = 1:size(LR_LF,4)
    for j = 1:size(LR_LF,5)
        if i == 5 && j == 5
            %-- do nothing
        else
            %--- Display the progress of the sift-flow computation
            msg = sprintf('  SIFTflow computation: %6.2f%%', k/(size(LR_LF,4)*size(LR_LF,5))*100);
            fprintf('%s',msg);
            lengthLastMsg = length(msg);
            pause(0.005);

            % Load the current sub-aperture image
            It = rgb2gray(LR_LF(:,:,:,i,j));
            
            % Extract the sift features for the target
            sift2 = mexDenseSIFT(It ,cellsize,gridspacing);
            
            % Compute the sift-flow
            [u_(:,:,i,j),v_(:,:,i,j)]=SIFTflowc2f(sift1,sift2,SIFTflowpara);
            
            % Correct the polarity of the flow vectors
            u_(:,:,i,j) = -u_(:,:,i,j);
            v_(:,:,i,j) = -v_(:,:,i,j);
            
            % Increment counter
            k = k + 1;
            
            %--- Clear the last entry
            fprintf(repmat('\b', 1, lengthLastMsg));
        end
    end
end