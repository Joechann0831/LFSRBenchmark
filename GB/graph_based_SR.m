function SR_LF = graph_based_SR(LR_LF,mf,sz,dispMax,sub_patch_size)
% This function is the interface function to use the graph based super
% resolution of Rossi et al.
%

%-- Configuration
factor = mf;
% Spatial resolution for each LR sub light field.
yResSubLR = sub_patch_size;
xResSubLR = sub_patch_size; 
% Pixel overlap between two adjacent LR sub light fields.
overlapLR = ceil(yResSubLR / 2);
% Merging window parameter.
alpha = 0.15;
% Set the PARFOR pool size.
poolSize = 4;

% Encode the light field into a cell array as used in their demo
ZLR = graph_based_encoding(LR_LF);

% Spatial resolution LR.
yResLR = size(ZLR{1, 1}, 1);
xResLR = size(ZLR{1, 1}, 2);
vRes = size(LR_LF,4);
hRes = size(LR_LF,5);

% Inter (views) graph parameters:
patRad = 3;             % Square patch radius.
intSigma = 0.7229;      % Weights decay.
% dispMax = 3;            % Maximum disparity.

% Energy function multipliers.
lambda0 = 1.0;
lambda1 = 0.2;
lambda2 = 0.0055;

% Number of inner and outer iterations. 
innerCycles = 200;
outerCycles = 2;

% Activation flag for the initial guess.
guessFlag = 1;

% Spatial resolution for each LR sub light field.
% yResSubLR = 100;
% xResSubLR = 100; 
warpMode = 'SQ';


% Super-resolve the light field ZLR. The computation is performed with pixels in double [0,1].
[ZInit, ZGB, ZInit_L, ZGB_L] = super( ...
    col2lf(double(lf2col(ZLR)) / 255, vRes, hRes, yResLR, xResLR, 3), ...
    factor, ...
    yResSubLR, xResSubLR, ...
    overlapLR, ...
    patRad, intSigma, dispMax, ...
    lambda0, lambda1, lambda2, ...
    innerCycles, outerCycles, ...
    warpMode, ...
    guessFlag, ...
    poolSize, ...
    alpha);

% Encode the light field into a cell array as used in their demo
SR_LF = graph_based_decoding(ZGB);

% Ensure that the super-resolved light field has the right dimensions
% SR_LF = SR_LF(1:sz(1),1:sz(2),:,:,:);


function X = graph_based_decoding(Z)

vRes = size(Z,1);
hRes = size(Z,2);

% Read the views paying attention to STANFORD reference system.
X = zeros(size(Z{1},1),size(Z{1},2),size(Z{1},3), vRes, hRes);
for t = 1:1:vRes
    for s = 1:1:hRes
        
        X(:,:,:,t,s) = Z{t, s};
    
    end
end
X = uint8(X*255);


function Z = graph_based_encoding(X)

vRes = size(X,4);
hRes = size(X,5);

% Read the views paying attention to STANFORD reference system.
Z = cell(vRes, hRes);
for t = 1:1:vRes
    for s = 1:1:hRes
        
        Z{t, s} = X(:,:,:,t,s);
    
    end
end
