function x=reshape(x,varargin)
[n,m]=size(x.values);
p=reshape(reshape(reshape(1:n*m,n,m),varargin{1},varargin{2}),n*m,1);
x.derivatives=x.derivatives(p,:);
x.values=reshape(x.values,varargin{1},varargin{2});