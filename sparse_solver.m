function [G1, impact] = sparse_solver(g1,psi,n_v,n_g)
%% NOTE: THIS SOLVER WILL BE UPDATED TO MATCH SCHUR_SOLVER EVENTUALLY
% Solves the rational expectation model with sparse iteration
%    method. You should use this method of the number of stable and
%    unstable solution differ a lot. Otherwise, full matrix 
%    decomposition method should be used.
%
% by SeHyoun Ahn, June 2016
%
% PARAMETERS/OUTPUTS:
%     Read attached documentation
%
% SYNTAX:
% [G1, impact] = sparse_solver(g1,psi,n_v,n_g)


%-100 is set to find negative eigenvalues. This value can be changed
[x,v,flag] = eigs(g1,n_g,-100);
impact = real(x*(x(n_v+1:end,:)\psi(n_v+1:end,:)));
G1 = real(x*v*(x(n_v+1:end,:)\[sparse(n_g,n_v),speye(n_g,n_g)]));
