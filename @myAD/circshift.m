function x = circshift(x,varargin)
% by SeHyoun Ahn, July 2016

aux = size(x.values);
locs = reshape(1:prod(aux),aux);
locs = circshift(locs,varargin{:});
x.derivatives = x.derivatives(locs(:),:);
x.values = circshift(x.values,varargin{:});