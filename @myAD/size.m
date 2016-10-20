function varargout = size(x, varargin)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
if nargin == 1
    [sx, sy] = size(x.values);
    if nargout <= 1
        varargout = {[sx, sy]};
    else
        varargout = {sx, sy};
    end
else
    sx = size(x.values, varargin{:});
    varargout = {sx};
end
