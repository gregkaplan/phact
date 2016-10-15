function y = subsasgn(y, S, x)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

if isa(x, 'myAD')
    if isa(y, 'myAD')
        y.values(S.subs{:}) = x.values;
        aux=size(y.values);
        locs=reshape(1:prod(aux),aux);
        locs=locs(S.subs{:});
        y.derivatives(locs(:),:) = x.derivatives;
    else
        aux=size(y);
        locs=reshape(1:prod(aux),aux);
        locs=locs(S.subs{:});
        l = size(x.derivatives,2);
        y = myAD(y,sparse(prod(aux),l));
        y.values(S.subs{:}) = x.values;
        y.derivatives(locs(:),:) = x.derivatives;
    end
else
    aux=size(y.values);
    locs=reshape(1:prod(aux),aux);
    locs=locs(S.subs{:});
    y.values(S.subs{:}) = x;
    y.derivatives(locs(:),:) = 0;
end