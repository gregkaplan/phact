function x = sin(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
x.derivatives = valXder(cos(x.values), x.derivatives);
x.values = sin(x.values);
