function [state_red,inv_state_red,n_g] = stateSpaceReduction(g0,g1,n_v,n_g,n_p,n_Z,varargin)
% Reduce Dimensionality of State Space given the direction to keep
%
% by SeHyoun Ahn, Sept 2016
%
% INPUTS: 
%   g0 = LHS matrix, only used to check it satisfies require form
%   g1 = Dynamics matrix
%   n_v = number of jump variables
%   n_g = number of state variables
%   n_p = number of static constraints
%   n_Z = number of exogenous variables
%
% OUTPUTS:
%   state_red = tranformation to get from full grid to reduced states
%   inv_state_red = inverse transform
%   n_g = number of state variables after reduction
%

%% Check to make sure that LHS satisfies the required form
if (nargin == 7)
  n_max_rank = varargin{1};
else
  n_max_rank = n_g;
end

loc = [1:n_v+n_g,n_v+n_g+n_p+1:n_v+n_g+n_p+n_Z];
if (max(max(abs(g0(loc,loc)-speye(n_v+n_g+n_Z))))~=0)
   error('Make sure that g0 is normalized.');
end

n_total = n_v + n_g;

%% Slice Dynamics Equation into Different Parts
B_pv = g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,1:n_v);
B_pZ = g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+n_p+1:n_v+n_g+n_p+n_Z);
B_pg = g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,n_v+1:n_v+n_g);
B_gg = g1(n_v+1:n_v+n_g,n_v+1:n_v+n_g);
B_gp = g1(n_v+1:n_v+n_g,n_v+n_g+1:n_v+n_g+n_p);

%% Orthonormalize B_pg as the direction to keep
[~,d0,V_g] = svd(full(B_pg),'econ');
aux = diag(d0);
n_Bpg = sum(aux>10*eps*aux(1));
V_g = V_g(:,1:n_Bpg);

%% Main loop to find State Space Reduction
n_old = 0;
n_new = n_Bpg;
V_current = V_g;
while (n_new>n_old) && (n_new < n_max_rank)
  n_old = n_new;
  V_current = B_gg'*V_current - B_pg'*(B_gp'*V_current);
  V_g = [V_g,full(V_current)];
  [~,d0,~] = svd(V_g,'econ');
  aux = diag(d0);
  n_new = sum(aux>10*eps*aux(1));
end
[V_g,d0,~] = svd(V_g,'econ');
aux = diag(d0);
n_g = sum(aux > 10*eps*aux(1));
V_g = V_g(:,1:n_g);

%% Build State-Space Reduction tranform
state_red = sparse(n_v+n_g+n_Z,n_total+n_p+n_Z);
state_red(1:n_v,1:n_v) = speye(n_v);
state_red(n_v+1:n_v+n_g,n_v+1:n_total) = V_g';
state_red(:,n_total+1:n_total+n_p) = 0;
state_red(n_v+n_g+1:n_v+n_g+n_Z,n_total+n_p+1:n_total+n_p+n_Z) = speye(n_Z);

%% Build inverse transform
inv_state_red = sparse(n_total+n_p+n_Z,n_v+n_g+n_Z);
inv_state_red(1:n_v,1:n_v) = speye(n_v);
inv_state_red(n_v+1:n_total,n_v+1:n_g+n_v) = V_g;
inv_state_red(n_total+1:n_total+n_p,n_v+1:n_v+n_g) = -B_pg*V_g;
inv_state_red(n_total+1:n_total+n_p,1:n_v) = -B_pv;
inv_state_red(n_total+1:n_total+n_p,n_v+n_g+1:n_v+n_g+n_Z) = -B_pZ;
inv_state_red(n_total+n_p+1:n_total+n_p+n_Z,n_v+n_g+1:n_v+n_g+n_Z) = speye(n_Z);
