function x = acos(x)
% Edited by SeHyoun Ahn, May 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
x.derivatives = valXder(-1./sqrt(1-x.values(:).^2), x.derivatives);
x.values = acos(x.values);
