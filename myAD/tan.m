function x = tan(x)
% Edited by SeHyoun Ahn, May 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
x.derivatives = valXder(1./cos(x.values(:)).^2, x.derivatives);
x.values = tan(x.values);
