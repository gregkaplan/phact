function [state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0_sparse(g0,g1,c,pi,psi)
% Solves out static constraints of the linear model by solving out last n_p
%     variables with 0 rows in g0.
%
% Input/Output/Rerences: It will be faster to read the codes below
%
% by SeHyoun, June 2106
%
% [state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0_sparse(g0,g1,c,pi,psi)

n = size(g0,1);
tmp = (max(abs(g0),[],2)==0);
redundant = find(tmp);
%n_red = length(redundant);
keep=find(1-tmp);
n_keep = length(keep);

inv_state_red = sparse(n,n_keep);
inv_state_red(keep,:) = speye(n_keep);
inv_state_red(redundant,:) = -g1(redundant,redundant)\g1(redundant,keep);

state_red = sparse(n_keep,n);
state_red(:,keep)=speye(n_keep);

g0=state_red*g0*inv_state_red;
g1=state_red*g1*inv_state_red;
g1=g0\g1;
psi=g0\state_red*psi;
pi=g0\state_red*pi;
c=g0\state_red*c;
