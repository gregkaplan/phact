function [state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0(g0,g1,c,pi,psi)
% Solves out static constraints of the linear model using SVD
%    Adapted from Chris Sims's gensys code
%
% Input/output/References: It will be faster to read the codes below
%
% by SeHyoun Ahn, June 2016
%
% SYNTAX:
% [state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0(g0,g1,c,pi,psi)
tmp=(max(abs([g0,psi]),[],2)==0);
redundant=find(tmp);
keep=find(1-tmp);
inv_state_red=sparse(null(full(g1(redundant,:))));

g0=inv_state_red'*g0*inv_state_red;
g1=inv_state_red'*g1*inv_state_red;
g1=g0\g1;
psi=g0\inv_state_red'*psi;
pi=g0\inv_state_red'*pi;
c=g0\inv_state_red'*c;

state_red = inv_state_red';
