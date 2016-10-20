function x = asin(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
x.derivatives = valXder(1./sqrt(1-x.values.^2), x.derivatives);
x.values = asin(x.values);
