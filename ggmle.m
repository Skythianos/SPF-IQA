function [mhat, ahat, bhat] = ggmle(x, options)
% GGMLE Parameter estimates for generalized Gaussian distributed data.
%
%       [MHAT, AHAT, BHAT] = GGMLE(X, OPTIONS) 
%
%       Returns the maximum likelihood estimates of the parameters of the 
%       generalized Gaussian distribution given the data in the vector X.
%
%	OPTIONS (option) is set by OPTIMSET to be used with FZERON
%
% Ref:  Minh N. Do and Martin Vetterli, "Wavelet-based texture retrieval using
%       generalized Gaussian density and Kullback-Leibler distance", IEEE Trans.
%       on Image Processing, Feb. 2002.
%
% Author: Minh N. Do, Dec. 1999


if min(size(x)) > 1
    error('The first argument in GGMLE must be a vector.');
end

if nargin < 2
    options = [];
end

x = x(:);

% Estimate the mean
mhat = mean(x);

% Initial estimate of beta by moment matching
dx = abs(x-mhat);
m1 = mean(dx);
m2 = mean(dx .* dx);
bhat = estbeta(m1, m2);

% Maximum likelihood estimates of beta and alpha
bhat = fzeron('dggbeta', bhat, options, dx);
ahat = (bhat * sum(dx .^ bhat) / length(x)) ^ (1 / bhat);
