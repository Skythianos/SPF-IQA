function r = ggrnd(m,a,b,n1,n2);
%GGRND Random matrices from generalized Gaussian distribution.
%   R = GGRND(M,A,B) returns a matrix of random numbers chosen   
%   from the generalized Gaussian distribution with parameters A and B.
%   The size of R is the common size of A and B if both are matrices.
%   If either parameter is a scalar, the size of R is the size of the other
%   parameter. Alternatively, R = GAMRND(M,A,B,N1,N2) returns an N1 by N2 matrix. 
%
%   PDF for generalized Gaussian distribution is defined as:
%   	F(X) = (B/(2*A*gamma(1/B))) * exp(-(abs(X-M)/A)^B)

%   GGRND uses the fact that if V is uniformly distributed on [-1, +1]
%   and Y is gamma(1+1/B,1) distributed then 
%                X = M + A*V*Y^(1/B) 
%   is generalized Gaussian distributed with mean M, scale A, and shape B
%
%   See also: GAMRND
%
% Author: Minh N. Do, Dec. 1999


%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986

if nargin < 3, 
   error('Requires at least three input arguments.'); 
end


if nargin == 3
   [errorcode rows columns] = rndcheck(3,3,m,a,b);
end

if nargin == 4
   [errorcode rows columns] = rndcheck(4,3,m,a,b,n1);
end

if nargin == 5
   [errorcode rows columns] = rndcheck(5,3,m,a,b,n1,n2);
end


if errorcode > 0
   error('Size information is inconsistent.');
end

% Initialize
lth = rows*columns;
m = m(:); a = a(:); b = b(:);

scalarm = (length(m) == 1);
if scalarm 
   m = m*ones(lth,1);
end

scalara = (length(a) == 1);
if scalara 
   a = a*ones(lth,1);
end

scalarb = (length(b) == 1);
if scalarb 
   b = b*ones(lth,1);
end

% V
v = 2 * (rand(lth,1) - 0.5);

% Y
y = gamrnd(1 + 1 ./ b, 1);

% R
r = m + a .* v .* y .^ (1 ./ b);

% Return NaN if a or b is not positive.
r(b <= 0 | a <= 0) = NaN;

r = reshape(r,rows,columns);