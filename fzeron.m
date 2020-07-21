function x = fzeron(fun, x0, varargin)
% FZERON Scalar nonlinear zero finding using Newton-Raphson method
%
%   X = FZERON(FUN,X0) tries to find a zero of FUN near X0. FUN (usually 
%   an M-file): FUN.M should take a scalar real value and return 2 real scalar
%   values when called with feval: [F,G]=feval(FUN,X); F is function
%   value and G is the derivative of the function evaluated at X.
%   The value X returned by FZERON is near a point where FUN changes
%   sign, or NaN if the search fails.
%
%   X = FZERON(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, MaxIter and TolX. Use OPTIONS = [] as a place
%   holder if no options are set.
%
%   X = FZERON(FUN,X0,OPTIONS,P1,P2,...) allows for additional arguments
%   which are passed to the function, [F,G]=feval(FUN,X,P1,P2,...).  
%   Pass an empty matrix for OPTIONS to use the default values.
%
%   See also: FZERO
%
% Author: Minh N. Do, Dec. 1999


% Initialization
% !!! OPTIMSET takes too much of time !!!
% defaultopt = optimset('display', 'none', 'TolX', 1e-6, 'MaxIter', 100);

if nargin < 2
    error('FZERON requires at least 2 inputs.');
elseif nargin == 2
    options = [];
else	% nargin >= 3
    options = varargin{1};
    varargin = varargin(2:end);
end
    
% options = optimset(defaultopt, options);
tol = optimget(options, 'TolX');
if isempty(tol)
    tol = 1e-6;
end

printtype = optimget(options, 'display');
if isempty(printtype)
    printtype = 'none';
end

maxiter = optimget(options, 'MaxIter');
if isempty(maxiter)
    maxiter = 100;
end

switch printtype
    case {'none', 'off'}
	trace = 0;
    case 'iter'
	trace = 2;
    case 'final'
	trace = 1;
    otherwise
	trace = 0;
end

if trace > 1 
    disp('   Func-count      x           f(x)');
end

if (~isfinite(x0) | length(x0) > 1)
    error('Second argument must be finite scalar.')
end

fcount = 1;
x = x0;
xsave = NaN;

while fcount <= maxiter
    [f, g] = feval(fun, x, varargin{:});    
    
    if trace > 1
	disp(sprintf('%5.0f   %13.6g %13.6g', fcount, x, f))
    end
    
    % Break for unnormal cases    
    if ~isfinite(f) | g == 0
	x = xsave;
	break;
    end
    
    % Save the last valid result
    if isfinite(f)
	xsave = x;
    end;
    
    % Iterative step
    delta = f / g;
    x = x - delta;

    if (abs(delta) <= tol)
	break;
    end

    fcount = fcount + 1;
end
