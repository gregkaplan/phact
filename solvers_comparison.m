% this is part of phact_solver

tmp=(max(abs([g0,psi]),[],2)==0); % identify the superfluous elements
redundant=find(tmp);
keep=find(1-tmp);
base=sparse(null(full(g1(redundant,:)))); % find the nullspace of the redundant parts

base = (base > 10^(-15)) .* base;

g0=base'*g0*base;
g1=base'*g1*base;
g1=g0\g1;
psi=g0\base'*psi;
pi=g0\base'*pi;
c=g0\base'*c;

%% VS

% this is cleanG0sparse

n = size(g0,1);

% tmp = (max(abs(g0),[],2)==0); % this line is an obvious candidate, but not the reason for the difference
tmp=(max(abs([g0,psi]),[],2)==0);

redundant = find(tmp);
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