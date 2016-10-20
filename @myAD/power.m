function x = power(x,y)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if isa(y, 'myAD')
    if isa(x, 'myAD')
        tmp1 = x.values.^(y.values);
        tmp2 = tmp1.*log(x.values);
        tmp3 = y.values.*x.values.^(y.values-1);
        x.derivatives = valXder(tmp3(:),x.derivatives) +...
            + valXder(tmp2(:),y.derivatives);
        x.values = tmp1;
    else
        y.values = x.^y.values;
        y.derivatives = valXder(y.values(:).*log(x(:)),y.derivatives);
        x = y;
    end
else
    x.derivatives = valXder(y(:).*x.values(:).^(y(:)-1),x.derivatives);
    x.values = x.values.^y;
end
