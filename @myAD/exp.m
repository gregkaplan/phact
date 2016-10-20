function x = exp(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
x.values = exp(x.values);
x.derivatives = valXder(x.values, x.derivatives);
