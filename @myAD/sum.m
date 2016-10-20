function x = sum(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
if nargin<2 || ~isa(varargin{1},'double')
    [n,m]=size(x.values);
    if n ==1
        x.values = sum(x.values,varargin{2:end});
        x.derivatives = sum(x.derivatives);
    else
        l=size(x.derivatives,2);
        x.derivatives = reshape(sum(reshape(x.derivatives,n,m*l)),m,l);
        x.values = sum(x.values,varargin{2:end});
    end
else
    x=sum(x',varargin{2:end})';
end
