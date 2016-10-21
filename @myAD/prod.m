function x = prod(x,varargin)
% by SeHyoun Ahn, July 2016

if nargin<2 || ~isa(varargin{1},'double')
    aux = size(x.values);
    reduction = find(aux>1);
    reduction = reduction(1);
    l=size(x.derivatives,2);
    
    aux1 = reshape(prod(reshape(kron(ones([ones(1,reduction-1),aux(reduction),ones(1,length(aux)-reduction)]), ...
        x.values),aux(reduction),prod(aux))),aux)./x.values;
    x.derivatives = reshape(sum(reshape(valXder(aux1,x.derivatives),aux(reduction),prod(aux)/aux(reduction)*l),varargin{:}),prod(aux)/aux(reduction),l);
    x.values = prod(x.values,varargin{:});
else
    aux = size(x.values);
    reduction = varargin{1};
    l=size(x.derivatives,2);
    p = permute(reshape(1:prod(aux),aux),[reduction,1:reduction-1,reduction+1:length(aux)]);
    aux1 = reshape(prod(reshape(kron(ones([aux(reduction),ones(1,length(aux)-1)]), ...
        x.values(p)),aux(reduction),prod(aux))),[aux(reduction),aux(1:reduction-1),aux(reduction+1:length(aux))])./x.values(p);
    x.derivatives = reshape(sum(reshape(valXder(aux1,x.derivatives(p,:)),aux(reduction),prod(aux)/aux(reduction)*l),varargin{2:end}),prod(aux)/aux(reduction),l);
    x.values = prod(x.values,varargin{:});
end