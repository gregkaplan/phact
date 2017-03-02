function [state_red,inv_state_red,n_g] = krylov_reduction(g0,g1,n_v,n_g,n_p,n_arnoldi,varargin)
% Reduce the dimensionality of state space using Krylov subspace method
%
% by SeHyoun Ahn, Sept 2016
%
% REFERENCES: <Our paper>
%
% PARAMETERS:
%    g0 = LHS matrix, only used to check it satisfies require form
%    g1 = Dynamics matrix
%    n_v = number of jump variables
%    n_g = number of state variables
%    n_p = number of static constraints
%    n_arnoldi = dimension of Krylov subspace
%    varargin(1) = other observables
%    varargin(2) = updating parameter for B_gv
%
% OUTPUTS:
%    state_red = tranformation to get from full grid to reduced states
%    inv_state_red = inverse transform
%    n_g = number of state variables after reduction


if nargin == 8
    observable = varargin{1};
    F = varargin{2};
elseif nargin == 7
    observable = varargin{1};
    F = [];
else
    observable = [];
    F = [];
end

%% Check to make sure that g0 is an identity matrix
%     This reducetion step assumes that g0 has been cleaned, so if
%     this part fails, make sure to run <clean_G0.m>
n_total = n_v + n_g;
if (max(max(abs(g0(1:n_total,1:n_total)-speye(n_total))))~=0)
    error('Make sure that g0 is normalized.');
end

% Slice Dynamics Equation into Different Parts
B_pv = g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,1:n_v);
B_pg = g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,n_v+1:n_v+n_g);
B_gg = g1(n_v+1:n_v+n_g,n_v+1:n_v+n_g);
B_gv = g1(n_v+1:n_v+n_g,1:n_v);
B_gp = g1(n_v+1:n_v+n_g,n_v+n_g+1:n_v+n_g+n_p);
% B_vp = g1(1:n_v,n_v+n_g+1:n_v+n_g+n_p);

%% Drop redundant Directions in B_pg
% In theory, since we do deflated block arnoldi, this step
%    might be redundant. Will be tested later.
obs = [B_pg;observable];
[~,d0,V_g] = svd(full(obs),'econ');
aux = diag(d0);
aux = aux/aux(1);
n_Bpg = sum(aux>10*eps);
V_g = bsxfun(@times,V_g(:,1:n_Bpg),aux(1:n_Bpg)');

% Compute Krylov Subspace
if isa(F,'function_handle')
    A = @(x) F(B_gv'*x)+B_gg'*x - B_pg'*(B_gp'*x);
elseif F
    A = @(x) F'*(B_gv'*x)+B_gg'*x - B_pg'*(B_gp'*x);
else
    A = @(x) B_gg'*x - B_pg'*(B_gp'*x);
end
[V_g,~] = deflated_block_arnoldi(A,V_g,n_arnoldi); 
n_g = size(V_g,2);

% Build State-Space Reduction transform
state_red = sparse(n_v+n_g,n_total+n_p);
state_red(1:n_v,1:n_v) = speye(n_v);
state_red(n_v+1:n_v+n_g,n_v+1:n_total) = V_g';
state_red(:,n_total+1:n_total+n_p) = 0;

% Build inverse transform
inv_state_red = sparse(n_total+n_p,n_v+n_g);
inv_state_red(1:n_v,1:n_v) = speye(n_v);
inv_state_red(n_v+1:n_total,n_v+1:n_g+n_v) = V_g;
inv_state_red(n_total+1:n_total+n_p,n_v+1:n_v+n_g) = -B_pg*V_g;
inv_state_red(n_total+1:n_total+n_p,1:n_v) = -B_pv;
