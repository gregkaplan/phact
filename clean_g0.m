function [base,inv_base,g1,c,pi,psi] = clean_g0(g0,g1,c,pi,psi)
% Inverts g0. For certain form on non-invertible g0 coming from constraints, it will compute the proper inverse.
% Written by SeHyoun, May 2106
n = size(g0,1);
tmp = (max(abs(g0),[],2)==0);
redundant = find(tmp);
%n_red = length(redundant);
keep=find(1-tmp);
n_keep = length(keep);

base = sparse(n,n_keep);
base(keep,:) = speye(n_keep);
base(redundant,:) = -g1(redundant,redundant)\g1(redundant,keep);

inv_base = sparse(n_keep,n);
inv_base(:,keep)=speye(n_keep);

g0=inv_base*g0*base;
g1=inv_base*g1*base;
g1=g0\g1;
psi=g0\inv_base*psi;
pi=g0\inv_base*pi;
c=g0\inv_base*c;

end
