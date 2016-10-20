function x=ctranspose(x)
% by SeHyoun Ahn, Jan 2016
% Only allows real derivatives, so this does not conjugate derivative
% values
[n,m]=size(x.values);
p=reshape(reshape(1:n*m,n,m)',n*m,1);
x.derivatives=x.derivatives(p,:);
x.values=x.values';
