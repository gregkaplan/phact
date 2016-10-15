function [g1,psi,pi,c,g0] = changeBasis(basis,inv_basis,g1,psi,pi,c,g0)
% Update dynamics with new basis
%
% by SeHyoun Ahn, June 2016

g1 = basis*g1*inv_basis;
pi = basis*pi;
psi = basis*psi;
c = basis*c;
if (nargin==7)
    g0= basis*g0*inv_basis;
else
    g0=speye(size(g1,1));
end
end
