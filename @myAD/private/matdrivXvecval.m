function output = matdrivXvecval(x,y)
[n,l]=size(x);
m=length(y);
n=n/m;
loc=reshape(1:n*m,n,m)';
output=reshape(sum(bsxfun(@times,y,reshape(x(loc(:),:),m,n*l))),n,l);