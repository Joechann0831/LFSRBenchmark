function LQ_LF = lf_awgn(HQ_LF,sig)
% This function is used to add additive white gaussion noise to the input
% light field. sig is a noise parameter.
%

% Initialize the low quality light field
LQ_LF = uint8(double(HQ_LF) +  sig*randn(size(HQ_LF)));