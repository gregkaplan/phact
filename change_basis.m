function [g1,psi,pi,c,g0] = change_basis(basis,inv_basis,g1,psi,pi,c,g0)
% Applies change of basis
%
% by SeHyoun Ahn, June 2016
%
% This is self-explanatory, so just check the codes below.
%
% Parameters: g0 is optional
%
% SYNTAX:
% [g1,psi,pi,c,g0] = change_basis(basis,inv_basis,g1,psi,pi,c,g0)


g1 = basis*g1*inv_basis;
pi = basis*pi;
psi = basis*psi;
c = basis*c;
if (nargin==7)
    g0= basis*g0*inv_basis;
else
    g0=speye(size(g1,1));
end
