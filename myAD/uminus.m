function x = uminus(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
x.values = - x.values;
x.derivatives = - x.derivatives;
