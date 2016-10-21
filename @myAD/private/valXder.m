function der = valXder(val, der)
%by SeHyoun Ahn, Jan 2016
der = bsxfun(@times, val(:),der);
