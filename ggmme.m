function [mu, alpha, beta] = ggmme(x)
% GGMME Moment matching estimator for generalized Gaussian pdf
%
%       [mu, alpha, beta] = ggmme(x)
%
% Input:
%       x: Input data
%
% Output:
%	mu, alpha, beta: paramaters of generalized Gaussian pdf
%	   f(x) = K * exp(-(abs(x-mu)/alpha)^beta)
%
% Author: Minh N. Do, Dec. 1999

x = x(:);

% Estimate the mean
mu = mean(x);

m1 = mean(abs(x - mu));
m2 = mean((x - mu) .* (x - mu));

% beta = F^(-1)(m1^2 / m2);
% alpha = m1 * gamma(1/beta) / gamma(2/beta);

fbeta = inline('exp(2 * gammaln(2 ./ x) - gammaln(3 ./ x) - gammaln(1 ./ x))', 'x');
F = sprintf('fbeta(x) - %g', m1^2 / m2);

try
    beta = fzero(F, [.01, 5]);
catch
    warning('m1^2 / m2 is out of the expected range');
    if (m1^2 / m2) > fbeta(5)
	beta = 5;
    else
	beta = 0.01;
    end    
end

alpha = m1 * exp(gammaln(1/beta) - gammaln(2/beta));
