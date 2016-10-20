function x = diff(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if nargin>1 && varargin{1}==2
    x = ctranspose(diff(ctranspose(x)));
else
    [n,m] = size(x.values);
    l = size(x.derivatives,2);
    x.values = diff(x.values);
    x.derivatives = reshape(x.derivatives,n,m*l);
    x.derivatives = diff(x.derivatives);
    x.derivatives = reshape(x.derivatives,(n-1)*m,l);
end
