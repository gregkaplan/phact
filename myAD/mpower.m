function x = mpower(x,y)
% by SeHyoun Ahn, Jan 2016
if numel(x)==1
    x=x.^y;
else
    error('Matrix power is not implemented. If elementwise exponentiation was intended use .^ ');
end
