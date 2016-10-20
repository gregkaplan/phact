function x = subsref(x, S)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
[n,m]=size(x.values);
locs=reshape(1:n*m,n,m);
locs=locs(S.subs{:});
x.derivatives = x.derivatives(locs(:),:);
x.values = x.values(S.subs{:});
