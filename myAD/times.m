function x = times(x,y)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
if isa(x, 'myAD')
    if isa(y, 'myAD')
        if max(size(y.values))==1
            x.derivatives = y.values.*x.derivatives + valXder(x.values(:),y.derivatives);
        elseif max(size(x.values))==1
            x.derivatives = valXder(y.values(:),x.derivatives) + x.values.*y.derivatives;
        else
            x.derivatives = valXder(y.values(:),x.derivatives) + valXder(x.values(:),y.derivatives);
        end
        x.values = x.values.*y.values;
    else
        if max(size(y))==1
            x.derivatives = y*x.derivatives;
        elseif max(size(x.values))==1
            x.derivatives = valXder(y(:), x.derivatives);
        else
            x.derivatives = valXder(y(:),x.derivatives);
        end
        x.values = x.values.*y;
    end
else
    if max(size(y.values))==1
        y.derivatives = valXder(x(:),y.derivatives);
    elseif max(size(x))==1
        y.derivatives = x.*y.derivatives;
    else
        y.derivatives = valXder(x(:),y.derivatives);
    end
    y.values = x.*y.values;
    x = y;
end
