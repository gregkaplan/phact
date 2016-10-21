function varargout = sort(x,varargin)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

if nargin==1 || varargin{1}==1
    aux = size(x.values);
    [x.values, idx] = sort(x.values,varargin{:});
    locs=repmat(reshape(aux(1)*(0:(prod(aux)/aux(1)-1)),[1,aux(2:end)]),[aux(1),ones(1,length(aux)-1)])+idx;
    x.derivatives = x.derivatives(locs(:),:);
else
    reduction = varargin{1};
    aux = size(x.values);
    p = permute(reshape(1:prod(aux),aux),[reduction,1:reduction-1,reduction+1:length(aux)]);
    S.type = '()';
    S.subs = {p};
    x = subsref(x,S);
    [x,idx] = sort(x,1,varargin{2:end});
    p = permute(reshape(1:prod(aux),[aux(reduction),aux(1:reduction-1),aux(reduction+1:length(aux))]),[2:reduction,1,reduction+1:length(aux)]);
    S.subs = {p};
    x = subsref(x,S);
    idx = idx(p);
end
varargout{1} = x;
varargout{2} = idx;