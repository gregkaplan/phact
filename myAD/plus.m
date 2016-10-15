function x = plus(x, y)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if isa(x, 'myAD')
    if isa(y, 'myAD')
        x.values = x.values + y.values;
        x.derivatives = x.derivatives + y.derivatives;
    else
        x.values = x.values + y;
    end
else
    y.values = x + y.values;
    x = y;
end
