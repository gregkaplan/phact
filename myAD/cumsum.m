function x = cumsum(x,varargin)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
aux = size(x.values);
l = size(x.derivatives,2);
if nargin == 1
    reduction = find(aux>1);
    reduction = reduction(1);
    x.derivatives = reshape(x.derivatives,aux(reduction),prod(aux)/aux(reduction)*l);
    x.derivatives = cumsum(x.derivatives);
    x.derivatives = reshape(x.derivatives,prod(aux),l);
    x.values = cumsum(x.values);
else
    reduction = varargin{1};
    x.values = cumsum(x.values,varargin{:});
    p = permute(reshape(1:prod(aux),aux),[reduction,1:reduction-1,reduction+1:length(aux)]);
    x.derivatives = reshape(cumsum(reshape(x.derivatives(p,:),aux(reduction),prod(aux)/aux(reduction)*l),1,varargin{2:end}),prod(aux),l);
    p = permute(reshape(1:prod(aux),aux),[2:reduction,1,reduction+1:length(aux)]);
    x.derivatives = x.derivatives(p,:);
end