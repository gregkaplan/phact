function x = diff(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
aux = size(x.values);
l = size(x.derivatives,2);
if nargin == 1
    reduction = find(aux>1);
    reduction = reduction(1);
    x.values = diff(x.values);
    x.derivatives = reshape(x.derivatives,aux(reduction),prod(aux)/aux(reduction)*l);
    x.derivatives = diff(x.derivatives);
    x.derivatives = reshape(x.derivatives,(aux(reduction)-1)*prod(aux)/aux(reduction),l);
elseif nargin == 2
    reduction = find(aux>1);
    reduction = reduction(1);
    x.values = diff(x.values);
    x.derivatives = reshape(x.derivatives,aux(reduction),prod(aux)/aux(reduction)*l);
    x.derivatives = diff(x.derivatives);
    x.derivatives = reshape(x.derivatives,(aux(reduction)-1)*prod(aux)/aux(reduction),l);
    if varargin{1} > 1
        x = diff(x,varargin{1}-1);
    end
else
    reduction = varargin{2};
    p = permute(reshape(1:prod(aux),aux),[reduction,1:reduction-1,reduction+1:length(aux)]);
    x.derivatives = reshape(x.derivatives(p,:),aux(reduction),prod(aux)/aux(reduction)*l);
    x.derivatives = diff(x.derivatives);
    x.derivatives = reshape(x.derivatives,(aux(reduction)-1)*prod(aux)/aux(reduction),l);
    p = permute(reshape(1:(prod(aux)/aux(reduction)*(aux(reduction)-1)),[aux(reduction)-1,aux(1:reduction-1),aux(reduction+1:length(aux))]),[2:reduction,1,reduction+1:length(aux)]);
    x.derivatives = x.derivatives(p,:);
    x.values = diff(x.values,[],varargin{2});
    if isa(varargin{1},'double') && varargin{1} >1
        x = diff(x,varargin{1}-1,varargin{2});
    end
end