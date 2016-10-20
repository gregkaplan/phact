function varargout = get(x, varargin)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if nargin <2
    varargout = {};
    display(x, inputname(1));
    return;
end

try
    varargout{1} = x.(varargin{1});
catch
    error(sprintf('Field %s doesn''t exist.', varargin{1}));
end
