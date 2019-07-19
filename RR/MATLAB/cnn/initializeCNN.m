function net = initializeCNN(convd,data)

%--------------------------------------------------------------------------
% Setup
%--------------------------------------------------------------------------
% setup the matconvnet library
vl_setupnn();

% Configuration
opts.useGpu = false ;                   % GPU enabled or not
opts.verbose = false ;                  % verbose or not

% Test that the convolution function works properly
try
  vl_nnconv(single(1),single(1),[]) ;
catch
  warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
  vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
               'enableImreadJpeg', false) ;
end
% Test that the GPU is working properly (if enabled)
if opts.useGpu
  try
    vl_nnconv(gpuArray(single(1)),gpuArray(single(1)),[]) ;
  catch
    warning('GPU support does not seem to be compiled in MatConvNet. Trying to compile it now.') ;
    vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
                 'enableImreadJpeg', false) ;
  end
end

% Initialize random bumber (ensure backward compatability)
if verLessThan('matlab','7.12')
  % MATLAB R2010b did not have rng()
  randn('state',0) ;
  rand('state',0) ;
else
  rng(0)  ;
end

% The EC2 has incorrect screen size which
% leads to a tiny font in figures

[~, hostname] = system('hostname') ;
if strcmp(hostname(1:3), 'ip-')
  set(0, 'DefaultAxesFontSize', 30) ;
end

%--------------------------------------------------------------------------
% Initialize CNN
%--------------------------------------------------------------------------
% Determine the dimension of each data element
M = data(1); N = data(2); C = data(3);

% Initialize the input layer
net.meta.inputSize = [M N C 1] ;

net.layers = { } ;

% Determine the number of convolutional layers
Nlayers = size(convd,2);

% We are using 9 neighbours -> 27 channels
nchannels = C;
for n = 1:Nlayers-1
    % Determine the filter size
    [filter_size, n_filters] = parse_conv(convd{n});
    net.layers{end+1} = struct(...
        'name', sprintf('conv%d',n), ...
        'type', 'conv', ...
        'weights', {xavier(filter_size,filter_size,nchannels,n_filters)}, ...
        'pad', 0, ...
        'stride', 1, ...
        'learningRate', [1 1], ...
        'weightDecay', [1 0]) ;

    net.layers{end+1} = struct(...
        'name', sprintf('relu%d',n), ...
        'type', 'relu') ;
    % Update the number of channels to the number of filters
    nchannels = n_filters;
end

% Determine the filter size
[filter_size, ~] = parse_conv(convd{end});

net.layers{end+1} = struct(...
  'name', 'prediction', ...
  'type', 'conv', ...
  'weights', {xavier(filter_size,filter_size,nchannels,C)}, ...
  'pad', 0, ... %(filter_size-1)/2, ...
  'stride', 1, ...
  'learningRate', [1 1], ...
  'weightDecay', [1 0]) ;

% net.layers{end+1} = struct(...
%   'name', 'conv1', ...
%   'type', 'conv', ...
%   'weights', {xavier(3,3,3,32)}, ...
%   'pad', 1, ...
%   'stride', 1, ...
%   'learningRate', [1 1], ...
%   'weightDecay', [1 0]) ;
% 
% net.layers{end+1} = struct(...
%   'name', 'relu1', ...
%   'type', 'relu') ;
% 
% net.layers{end+1} = struct(...
%   'name', 'conv2', ...
%   'type', 'conv', ...
%   'weights', {xavier(3,3,32,32)}, ...
%   'pad', 1, ...
%   'stride', 1, ...
%   'learningRate', [1 1], ...
%   'weightDecay', [1 0]) ;
% 
% net.layers{end+1} = struct(...
%   'name', 'relu2', ...
%   'type', 'relu') ;

% net.layers{end+1} = struct(...
%   'name', 'prediction', ...
%   'type', 'conv', ...
%   'weights', {xavier(3,3,32,3)}, ...
%   'pad', 1, ...
%   'stride', 1, ...
%   'learningRate', [1 1], ...
%   'weightDecay', [1 0]) ;

% Consolidate the network, fixing any missing option
% in the specification above

net = vl_simplenn_tidy(net) ;

function weights = xavier(varargin)
%XAVIER  Xavier filter initialization.
%   WEIGHTS = XAVIER(H, W, C, N) initializes N filters of support H x
%   W and C channels using Xavier method. WEIGHTS = {FILTERS,BIASES}is
%   a cell array containing both filters and biases.
%
% See also:
% Glorot, Xavier, and Yoshua Bengio.
% "Understanding the difficulty of training deep feedforward neural networks."
% International conference on artificial intelligence and statistics. 2010.

filterSize = [varargin{:}] ;
scale = sqrt(2/prod(filterSize(1:3))) ;
filters = randn(filterSize, 'single') * scale ;
biases = zeros(filterSize(4),1,'single') ;
weights = {filters, biases} ;


function [filter_size, n_filters] = parse_conv(convd)

idx1 = strfind(convd,'conv') + length('conv');
idx2  = strfind(convd,'-');

filter_size = str2double(convd(idx1:idx2-1));
n_filters   = str2double(convd(idx2+1:end));

