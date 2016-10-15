function x=repmat(x,varargin)
[n,m]=size(x.values);
p=reshape(repmat(reshape(1:n*m,n,m),varargin{1},varargin{2}),n*m*varargin{1}*varargin{2},1);
x.derivatives=x.derivatives(p,:);
x.values=repmat(x.values,varargin{1},varargin{2});