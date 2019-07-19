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

opts.expDir = fullfile('data','exp') ;
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
