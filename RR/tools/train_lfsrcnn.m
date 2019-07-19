function train_lfsrcnn(mf)
% This function is used to train the lfsrcnn light field super-resolution
% models

clc; close all;

%---- VLFeat library
addpath('../MATLAB/matconvnet/matlab/');           % MATCONVNET scripts
addpath('../MATLAB/matconvnet/matlab/simplenn/');  % Simple NN scripts
addpath('../MATLAB/cnn/');
%addpath('C:\Users\rrfarr\Documents\LR-CNN-LFSR\MATLAB\matconvnet\examples\mnist');
% Seed the random number generator so that we can repeat it
rng('default');

trainOpts.batchSize         = 128;      % batch size
trainOpts.learningRate      = 0.0001;   % learning rate
trainOpts.plotDiagnostics   = false;    % set default value
trainOpts.numEpochs         = 2000;       % number of epochs
trainOpts.errorFunction     = 'none' ;  % set default value
trainOpts.gpus = 1;
transferModelFilename = [];
patch_size = 64;

% Define the neural network model
conv_data = {'conv9-64','conv1-32','conv5-1'};

% Derive the path where the network model will be stored at each epoch    
trainOpts.expDir = sprintf('../DATA/superresolution/lfsrcnn/model/mf-%d/',mf); % Specify path where data

% Define where the imdb data will be stored
imdb_filename  = sprintf('../DATA/superresolution/lfsrcnn/imdb_x%d.mat',mf);

% Initialize the CNN
net = initializeCNN(conv_data,[patch_size,patch_size,1]);

% Add a loss (using a custom layer)
net = addLossLayer(net, @l2LossForward, @l2LossBackward);

% Load the imdb data
load(imdb_filename);

if ~exist(trainOpts.expDir,'dir')
    mkdir(trainOpts.expDir);
end

% Train the cnn based on the input configuration
cnn_train(net, imdb, @getBatch, transferModelFilename, trainOpts);

function net = initializeCNN(convd,data)

%--------------------------------------------------------------------------
% Setup
%--------------------------------------------------------------------------
% setup the matconvnet library
vl_setupnn();

% Configuration
opts.useGpu = false ;                   % GPU enabled or not
opts.verbose = false ;                  % verbose or not

try
    vl_nnconv(gpuArray(1),gpuArray(1),[]) ;
catch
    vl_compilenn('enableGpu', true, ...
        'cudaRoot', 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0', ...
        'cudaMethod', 'nvcc', ...
        'enableCudnn', true, ...
        'cudnnRoot', 'C:\cudnn-7.0-win-x64-v4.0\cuda');
end

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

cnnmodel = load('../DATA/superresolution/lfsrcnn/lf-srcnn-init.mat');

net.layers = cnnmodel.net.layers;

% Determine the number of convolutional layers
Nlayers = size(convd,2);

% We are using 9 neighbours -> 27 channels
nchannels = C;
% for n = 1:Nlayers-1
%     % Determine the filter size
%     [filter_size, n_filters] = parse_conv(convd{n});
%     
%     %a = model.net.layers{n}.weights{1};
%     net.layers{end+1} = struct(...
%         'name', sprintf('conv%d',n), ...
%         'type', 'conv', ...
%         'weights', cnnmodel.net.layers{2*(n-1)+1}.weights, ...
%         'pad', 0, ...
%         'stride', 1, ...
%         'learningRate', [1 1], ...
%         'weightDecay', [1 0]) ;
% 
%     net.layers{end+1} = struct(...
%         'name', sprintf('relu%d',n), ...
%         'type', 'relu') ;
%     % Update the number of channels to the number of filters
%     nchannels = n_filters;
% end
% 
% % Determine the filter size
% [filter_size, ~] = parse_conv(convd{end});
% 
% net.layers{end+1} = struct(...
%   'name', 'prediction', ...
%   'type', 'conv', ...
%   'weights', cnnmodel.net.layers{2*(n-1)+1}.weights, ...
%   'pad', 0, ... %(filter_size-1)/2, ...
%   'stride', 1, ...
%   'learningRate', [1 1], ...
%   'weightDecay', [1 0]) ;

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


function [net, info] = cnn_train(net, imdb, getBatch, transferModelFilename, varargin)

%CNN_TRAIN  An example implementation of SGD for training CNNs
%    CNN_TRAIN() is an example learner implementing stochastic
%    gradient descent with momentum to train a CNN. It can be used
%    with different datasets and tasks by providing a suitable
%    getBatch function.
%
%    The function automatically restarts after each training epoch by
%    checkpointing.
%
%    The function supports training on CPU or on one or more GPUs
%    (specify the list of GPU IDs in the `gpus` option). Multi-GPU
%    support is relatively primitive but sufficient to obtain a
%    noticable speedup.

% Copyright (C) 2014-15 Andrea Vedaldi.
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).

opts.expDir = varargin{1}.expDir; 
opts.continue = true ;
opts.batchSize = 256 ;
opts.numSubBatches = 1 ;
opts.train = [] ;
opts.val = [] ;
opts.gpus = [] ;
opts.prefetch = false ;
opts.numEpochs = 20 ;
opts.learningRate = 0.001 ;
opts.weightDecay = 0.0005 ;
opts.momentum = 0.9 ;
opts.memoryMapFile = fullfile(tempdir, 'matconvnet.bin') ;
opts.profile = false ;

opts.conserveMemory = true ;
opts.backPropDepth = +inf ;
opts.sync = false ;
opts.cudnn = true ;
opts.errorFunction = 'multiclass' ;
opts.errorLabels = {} ;
opts.plotDiagnostics = false ;
opts.plotStatistics = true;
opts = vl_argparse(opts, varargin) ;

if ~exist(opts.expDir, 'dir'), mkdir(opts.expDir) ; end
if isempty(opts.train), opts.train = find(imdb.images.set==1) ; end
if isempty(opts.val), opts.val = find(imdb.images.set==2) ; end
if isnan(opts.train), opts.train = [] ; end
if isnan(opts.val), opts.val = [] ; end

% -------------------------------------------------------------------------
%                                                    Network initialization
% -------------------------------------------------------------------------

net = vl_simplenn_tidy(net); % fill in some eventually missing values
net.layers{end-1}.precious = 1; % do not remove predictions, used for error
vl_simplenn_display(net, 'batchSize', opts.batchSize) ;

evaluateMode = isempty(opts.train) ;

%--- Initialize the network parameters
if ~evaluateMode
    for i=1:numel(net.layers) % Loop across all layers
        if isfield(net.layers{i}, 'weights')
            % Determine the number of weight components
            J = numel(net.layers{i}.weights) ;
            % Initialize layer i and the momentum j to all zeros
            for j=1:J
                net.layers{i}.momentum{j} = zeros(size(net.layers{i}.weights{j}), 'single') ;
            end
            % Set the learning rate for each layer to 1
            if ~isfield(net.layers{i}, 'learningRate')
                net.layers{i}.learningRate = ones(1, J, 'single') ;
            end
            % Set the weight decay parameter to 1
            if ~isfield(net.layers{i}, 'weightDecay')
                net.layers{i}.weightDecay = ones(1, J, 'single') ;
            end
        end
    end
end

%--- Setup the GPUs
numGpus = numel(opts.gpus) ;
if numGpus > 1
    if isempty(gcp('nocreate')),
        parpool('local',numGpus) ;
        spmd, gpuDevice(opts.gpus(labindex)), end
    end
elseif numGpus == 1
    gpuDevice(opts.gpus)
end
if exist(opts.memoryMapFile), delete(opts.memoryMapFile) ; end

%--- Setup error calculation function
hasError = true ;
if ischar(opts.errorFunction)
    switch opts.errorFunction
        case 'none'
        opts.errorFunction = @error_none ;
        hasError = false ;
    case 'multiclass'
        opts.errorFunction = @error_multiclass ;
        if isempty(opts.errorLabels), opts.errorLabels = {'top1err', 'top5err'} ; end
    case 'binary'
        opts.errorFunction = @error_binary ;
        if isempty(opts.errorLabels), opts.errorLabels = {'binerr'} ; end
    otherwise
        error('Unknown error function ''%s''.', opts.errorFunction) ;
    end
end

% -------------------------------------------------------------------------
%                                                        Train and validate
% -------------------------------------------------------------------------
% Specify the filename of the model per epoch to be stored
modelPath = @(ep) fullfile(opts.expDir, sprintf('net-epoch-%d.mat', ep));
% Define the training figure
modelFigPath = fullfile(opts.expDir, 'net-train.pdf') ;


start = opts.continue * findLastCheckpoint(opts.expDir) ;
if start >= 1
    fprintf('%s: resuming by loading epoch %d\n', mfilename, start) ;
    load(modelPath(start), 'net', 'info') ;
    net = vl_simplenn_tidy(net) ; % just in case MatConvNet was updated
elseif ~isempty(transferModelFilename)
    load(transferModelFilename);
    % Clear the info structure
    clearvars info;
    net = vl_simplenn_tidy(net) ;
    fprintf('%s: Using a pre-trained model\n', mfilename) ;
end

for epoch=start+1:opts.numEpochs

  % train one epoch and validate
  learningRate = opts.learningRate(min(epoch, numel(opts.learningRate))) ;
  % Training data is randomized
  train = opts.train(randperm(numel(opts.train))) ; % shuffle
  % Consider all testing data
  val = opts.val ;

  if numGpus <= 1
      % Process an epoch on the training data
      [net,stats.train,prof] = process_epoch(opts, getBatch, epoch, train, ...
          learningRate, imdb, net) ;

      [~,stats.val] = process_epoch(opts, getBatch, epoch, val, 0, imdb, net) ;
      if opts.profile
          profile('viewer') ;
          keyboard ;
      end
  else
        fprintf('%s: sending model to %d GPUs\n', mfilename, numGpus) ;
        spmd(numGpus)
        [net_, stats_train_,prof_] = process_epoch(opts, getBatch, epoch, train, learningRate, imdb, net) ;
        [~, stats_val_] = process_epoch(opts, getBatch, epoch, val, 0, imdb, net_) ;
    
        end
        net = net_{1} ;
        stats.train = sum([stats_train_{:}],2) ;
        stats.val = sum([stats_val_{:}],2) ;
        if opts.profile
            mpiprofile('viewer', [prof_{:,1}]) ;
            keyboard ;
        end
        clear net_ stats_train_ stats_val_ ;
  end
  
  % save
  if evaluateMode, sets = {'val'} ; else sets = {'train', 'val'} ; end
    for f = sets
        f = char(f) ;
        n = numel(eval(f)) ;
        info.(f).speed(epoch) = n / stats.(f)(1) * max(1, numGpus) ;
        info.(f).objective(epoch) = stats.(f)(2) / n ;
        info.(f).error(:,epoch) = stats.(f)(3:end) / n ;
    end
    if ~evaluateMode
        fprintf('%s: saving model for epoch %d\n', mfilename, epoch) ;
        tic ;
        save(modelPath(epoch), 'net', 'info') ;
        fprintf('%s: model saved in %.2g s\n', mfilename, toc) ;
    end

    if opts.plotStatistics
        switchfigure(1) ; clf ;
        subplot(1,1+hasError,1) ;
        if ~evaluateMode
            semilogy(1:epoch, info.train.objective, '.-', 'linewidth', 2) ;
            hold on ;
        end
        semilogy(1:epoch, info.val.objective, '.--') ;
        xlabel('training epoch') ; ylabel('energy') ;
        grid on ;
        h=legend(sets) ;
        set(h,'color','none');
        title('objective') ;
        if hasError
            subplot(1,2,2) ; leg = {} ;
            if ~evaluateMode
                plot(1:epoch, info.train.error', '.-', 'linewidth', 2) ;
                hold on ;
                leg = horzcat(leg, strcat('train ', opts.errorLabels)) ;
            end
            plot(1:epoch, info.val.error', '.--') ;
            leg = horzcat(leg, strcat('val ', opts.errorLabels)) ;
            set(legend(leg{:}),'color','none') ;
            grid on ;
            xlabel('training epoch') ; ylabel('error') ;
            title('error') ;
        end
        drawnow ;
        print(1, modelFigPath, '-dpdf') ;
    end
end

% -------------------------------------------------------------------------
function  [net_cpu,stats,prof] = process_epoch(opts, getBatch, epoch, subset, learningRate, imdb, net_cpu)
% -------------------------------------------------------------------------

% move the CNN to GPU (if needed)
numGpus = numel(opts.gpus) ;
if numGpus >= 1
  net = vl_simplenn_move(net_cpu, 'gpu') ;
  one = gpuArray(single(1)) ;
else
  net = net_cpu ;
  net_cpu = [] ;
  one = single(1) ;
end

% assume validation mode if the learning rate is zero
training = learningRate > 0 ;
if training
  mode = 'train' ;
  evalMode = 'normal' ;
else
  mode = 'val' ;
  evalMode = 'test' ;
end

% turn on the profiler (if needed)
if opts.profile
  if numGpus <= 1
    prof = profile('info') ;
    profile clear ;
    profile on ;
  else
    prof = mpiprofile('info') ;
    mpiprofile reset ;
    mpiprofile on ;
  end
end

res = [] ;
mmap = [] ;
stats = [] ;
start = tic ;

for t=1:opts.batchSize:numel(subset)
  fprintf('%s: epoch %02d: %3d/%3d: ', mode, epoch, ...
          fix((t-1)/opts.batchSize)+1, ceil(numel(subset)/opts.batchSize)) ;
  batchSize = min(opts.batchSize, numel(subset) - t + 1) ;
  numDone = 0 ;
  error = [] ;
  for s=1:opts.numSubBatches
    % get this image batch and prefetch the next
    batchStart = t + (labindex-1) + (s-1) * numlabs ;
    batchEnd = min(t+opts.batchSize-1, numel(subset)) ;
    batch = subset(batchStart : opts.numSubBatches * numlabs : batchEnd) ;
    [im, labels] = getBatch(imdb, batch) ;
    % Normalize the batch
    im = single(im)/255;
    labels = single(labels)/255;
    % Consider only the center part of the image because of conv the
    % dimention of the activation is reduced
    labels = labels(9:end-8,9:end-8,:,:);

    if opts.prefetch
      if s==opts.numSubBatches
        batchStart = t + (labindex-1) + opts.batchSize ;
        batchEnd = min(t+2*opts.batchSize-1, numel(subset)) ;
      else
        batchStart = batchStart + numlabs ;
      end
      nextBatch = subset(batchStart : opts.numSubBatches * numlabs : batchEnd) ;
      getBatch(imdb, nextBatch) ;
    end

    if numGpus >= 1
      im = gpuArray(im) ;
    end

    % evaluate the CNN
    net.layers{end}.class = labels ;
    if training, dzdy = one; else, dzdy = [] ; end
    % Compute the feedforward part of the network
    res = vl_simplenn(net, im, dzdy, res, ...
                      'accumulate', s ~= 1, ...
                      'mode', evalMode, ...
                      'conserveMemory', opts.conserveMemory, ...
                      'backPropDepth', opts.backPropDepth, ...
                      'sync', opts.sync, ...
                      'cudnn', opts.cudnn) ;

    % accumulate training errors
    error = sum([error, [...
      sum(double(gather(res(end).x))) ;
      reshape(opts.errorFunction(opts, labels, res),[],1) ; ]],2) ;
    numDone = numDone + numel(batch) ;
  end % next sub-batch

  % gather and accumulate gradients across labs
  if training
    if numGpus <= 1
      [net,res] = accumulate_gradients(opts, learningRate, batchSize, net, res) ;
    else
      if isempty(mmap)
        mmap = map_gradients(opts.memoryMapFile, net, res, numGpus) ;
      end
      write_gradients(mmap, net, res) ;
      labBarrier() ;
      [net,res] = accumulate_gradients(opts, learningRate, batchSize, net, res, mmap) ;
    end
  end

  % collect and print learning statistics
  time = toc(start) ;
  stats = sum([stats,[0 ; error]],2); % works even when stats=[]
  stats(1) = time ;
  n = t + batchSize - 1 ; % number of images processed overall
  speed = n/time ;
  fprintf('%.1f Hz%s\n', speed) ;

  m = n / max(1,numlabs) ; % num images processed on this lab only
  fprintf(' obj:%.3g', stats(2)/m) ;
  for i=1:numel(opts.errorLabels)
    fprintf(' %s:%.3g', opts.errorLabels{i}, stats(i+2)/m) ;
  end
  fprintf(' [%d/%d]', numDone, batchSize);
  fprintf('\n') ;

  % collect diagnostic statistics
  if training & opts.plotDiagnostics
    switchfigure(2) ; clf ;
    diag = [res.stats] ;
    barh(horzcat(diag.variation)) ;
    set(gca,'TickLabelInterpreter', 'none', ...
      'YTickLabel',horzcat(diag.label), ...
      'YDir', 'reverse', ...
      'XScale', 'log', ...
      'XLim', [1e-5 1]) ;
    drawnow ;
  end

end

% switch off the profiler
if opts.profile
  if numGpus <= 1
    prof = profile('info') ;
    profile off ;
  else
    prof = mpiprofile('info');
    mpiprofile off ;
  end
else
  prof = [] ;
end

% bring the network back to CPU
if numGpus >= 1
  net_cpu = vl_simplenn_move(net, 'cpu') ;
else
  net_cpu = net ;
end
