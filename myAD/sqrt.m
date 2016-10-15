function x = sqrt(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
x.values = sqrt(x.values);
x.derivatives = valXder(0.5./x.values, x.derivatives);
