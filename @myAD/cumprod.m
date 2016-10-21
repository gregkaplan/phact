function x = cumprod(x,varargin)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if nargin ==3 && strcmp(varargin{2},'reverse')
    aux = size(x.values);
    reduction = varargin{1};
    p = sort(reshape(1:prod(aux),aux),varargin{1},'descend');
    S.type='()';
    S.subs={p};
    x = subsref(x,S);
    x = cumprod(x,reduction);
    x = subsref(x,S);
elseif nargin == 1 || (nargin>1 && varargin{1}==1)
    aux = size(x.values);
    reduction = find(aux>1);
    reduction = reduction(1);
    l=size(x.derivatives,2);
    S.type='()';
    derivatives = sparse(prod(aux),l);
    S.subs = {1,':'};
    derivatives(1+aux(reduction)*(0:(prod(aux)/aux(reduction)-1)),:) = getderivs(subsref(x,S));
    for j = 2:aux(reduction)
        S.subs={1:j,':'};
        derivatives(j+aux(reduction)*(0:(prod(aux)/aux(reduction)-1)),:) = getderivs(prod(subsref(x,S)));
    end
    x.derivatives = derivatives;
    x.values = cumprod(x.values);
else
    aux = size(x.values);
    reduction = varargin{1};
    x = permute(x,[reduction,1:(reduction-1),(reduction+1):length(aux)]);
    x = cumprod(x,1,varargin{2:end});
    x = permute(x,[2:reduction,1,reduction+1:length(aux)]);
end