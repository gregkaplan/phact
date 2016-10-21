function x = sqrt(x)
% Edited by SeHyoun Ahn, May 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
x.values = sqrt(x.values);
x.derivatives = valXder(0.5./x.values(:), x.derivatives);
