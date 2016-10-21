function x=spdiags(x,varargin)
% by SeHyoun Ahn, Jan 2016

% Note that spdiags will keep the derivative matrix in column-major order %
% Therefore, given matrix A(i,j), the first few rows of the derivative
% matrices are
% D_x A_{1,1} (These are row vectors of derivative of A_{1,1})
% D_x A_{2,1}
% ...
% D_x A_{n,1}
% D_x A_{1,2}
% ...
offset=varargin{1};
nrow = varargin{2};
ncol = varargin{3};
x.values=spdiags(x.values,offset,nrow,ncol);
l=size(x.derivatives,2);
if offset<0
    [i,j,v] = find(x.derivatives(1:min((nrow+offset),ncol),:));
    x.derivatives = sparse((1-offset)+(nrow+1)*(i-1),j,v,nrow*ncol,l);
else
    [i,j,v] = find(x.derivatives((offset+1):nrow,:));
    x.derivatives = sparse((offset*nrow)+1+(nrow+1)*(i-1),j,v,nrow*ncol,l);
end