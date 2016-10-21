function x=repmat(x,varargin)
% by SeHyoun Ahn, Jan 2016
aux=size(x.values);
p=repmat(reshape(1:prod(aux),aux),varargin{:});
x.derivatives=x.derivatives(p(:),:);
x.values=repmat(x.values,varargin{:});
