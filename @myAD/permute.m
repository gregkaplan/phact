function x=permute(x,order)
% by SeHyoun Ahn, July 2016
aux=size(x.values);
p=permute(reshape(1:prod(aux),aux),order);
x.derivatives=x.derivatives(p(:),:);
x.values=permute(x.values,order);
