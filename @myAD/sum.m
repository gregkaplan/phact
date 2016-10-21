function x = sum(x,varargin)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

if nargin<2 || ~isa(varargin{1},'double')
    aux = size(x.values);
    reduction = find(aux>1);
    reduction = reduction(1);
    l=size(x.derivatives,2);
    x.values = sum(x.values,varargin{:});
    x.derivatives = reshape(sum(reshape(x.derivatives,aux(reduction),prod(aux)/aux(reduction)*l),varargin{:}),prod(aux)/aux(reduction),l);
else
    aux = size(x.values);
    reduction = varargin{1};
    l=size(x.derivatives,2);
    x.values = sum(x.values,varargin{:});
    p = permute(reshape(1:prod(aux),aux),[reduction,1:reduction-1,reduction+1:length(aux)]);
    x.derivatives = reshape(sum(reshape(x.derivatives(p,:),aux(reduction),prod(aux)/aux(reduction)*l),varargin{2:end}),prod(aux)/aux(reduction),l);
end