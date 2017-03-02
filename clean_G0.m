function [base,inv_base,g1,c,pi,psi] = cleanG0(g0,g1,c,pi,psi)
% Inverts g0 using SVD.
%
% by SeHyoun Ahn, June 2016

tmp=(max(abs([g0,psi]),[],2)==0);
redundant=find(tmp);
keep=find(1-tmp);
base=sparse(null(full(g1(redundant,:))));

g0=base'*g0*base;
g1=base'*g1*base;
g1=g0\g1;
psi=g0\base'*psi;
pi=g0\base'*pi;
c=g0\base'*c;

inv_base = base';
