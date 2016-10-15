function x = sin(x)
% Edited by SeHyoun Ahn, May 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
x.derivatives = valXder(cos(x.values(:)), x.derivatives);
x.values = sin(x.values);
