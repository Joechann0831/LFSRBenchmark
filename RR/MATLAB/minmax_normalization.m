function [Y,min_X,max_X] = minmax_normalization(X)

min_X = repmat(min(X),[size(X,1),1]);
max_X = repmat(max(X),[size(X,1),1]);
Y = (X - min_X)./(max_X - min_X);
min_X = min(X); max_X = max(X);