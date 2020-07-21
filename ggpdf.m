function y = ggpdf(x, m, a, b)
% GGPDF	Generalized Gaussian probability density function
%	Y = GGPDF(X, M, A, B) returns the generalized Gaussian density function
%	with parameters A and B, at the values in X.
%
%	Y = (B/(2*A*gamma(1/B))) * exp(-(abs(X-M)/A)^B)
%
% Author: Minh N. Do, Dec. 1999

%   Return NaN if the arguments are outside their respective limits.
if (a <= 0 | b <= 0)
    tmp = NaN;
    y = tmp(ones(size(x)));    
else
    y = log(b) - log(2) - log(a) - gammaln(1/b) - (abs(x-m).^b) ./ (a^b);
    y = exp(y);
end