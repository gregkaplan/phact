function varargout = sort(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
if nargin>1 && isa(varargin{1},'double')
    if varargin{1}==1
        [n,m] = size(x.values);
        z = x;
        [z.values, idx] = sort(x.values);
        tmp = bsxfun(@plus,n*(0:m-1),idx);
        z.derivatives = x.derivatives(tmp,:);
    else
        [n,~] = size(x.values);
        z = x;
        [z.values, idx] = sort(x.values,varargin{:});
        tmp = bsxfun(@plus,(1:n)',n*(idx-1));
        z.derivatives = x.derivatives(tmp,:);
    end
else
    [n,m] = size(x.values);
    z = x;
    [z.values, idx] = sort(x.values,varargin{1:end});
    tmp = bsxfun(@plus,n*(0:m-1),idx);
    z.derivatives = x.derivatives(tmp,:);
end
varargout{1} = z;
varargout{2} = idx;
