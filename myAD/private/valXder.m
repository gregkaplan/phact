function der = valXder(val, der)
der = bsxfun(@times, val(:),der);
