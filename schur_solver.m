function [G1,C,impact,eu,F]=schur_solver(g0,g1,c,psi,pi,continuous,check_exist,check_uniq,varargin)
% Solver based on Schur decomposition (This is based on Sims's gensys, but adjusted
%    to take advantage of special case that applies to heterogeneous agent models).
%
% by SeHyoun Ahn, Dec 2016
%
% REFERENCES: Sims, Christopher A. "Solving linear rational expectations models."
%             Computational economics 20.1 (2002): 1-20.
%
% PARAMETERS:
%    Check Sims's paper for notation
%    continuous = 1 for a continuous time problem (default)
%                0 for a discrete time problem
%    check_exist = 1 check for existence of solution (default)
%                 0 existence has been checked before
%    check_uniq = 0 do not check for uniqueness (default)
%                1 check for uniqueness
%    varargin{1} = number of unstable roots 
%
% OUPUTS:
%    eu(1) = 1 solution exists
%            0 no stable solution
%           -5 Flag for checking existence is off
%    eu(2) = 1 solution is unique
%            0 solution is not unique
%           -5 Flag for checking uniqueness is off


%% Set Default Values %%
switch nargin
    case 9
        n_v = varargin{1};
    case 8
        n_v = -1;
    case 7
        check_uniq = 0;
        n_v = -1;
    case 6
        check_uniq = 0;
        check_exist = 1;
        n_v = -1;
    case 5
        check_uniq = 0;
        check_exist = 1;
        continuous = 1;
        n_v = -1;
    case 1:4
        error('Insufficient Number of Input Variables');
end

realsmall = sqrt(eps)*10;
eu = [-5;-5];
n = size(g1,1);

%% Schur Decomposition
[U,T] = schur(full(g1),'real');
if continuous
    [U,T] = ordschur(U,T,'lhp');
    g_eigs = ordeig(T);
    nunstab = sum(g_eigs>0);
    if n_v > -1 
        if nunstab > n_v
            aux = real(g_eigs);
            aux = sort(aux(aux>0),'descend');
            warning('There are more than n_v number of positive eigenvalues with smallest values:');
            disp(aux(end-nunstab+n_v+1:end));
        end
        nunstab = n_v;
    end
else
    [U,T] = ordschru(U,T,'udi');
    g_eigs = abs(ordeig(T));
    nunstab = sum(g_eigs<1);
    if n_v > -1
        if nunstab > n_v
            aux = abs(g_eigs);
            aux = sort(aux(aux>1),'descend');
            warning('There are more than n_v number of positive eigenvalues with smallest values:');
            disp(aux(end-nunstab+n_v+1:end));
        end
        nunstab = n_v;
    end
end

u2 = U(:,n-nunstab+1:n)';
u1 = U(:,1:n-nunstab)';

etawt = u2*pi;
[ueta,deta,veta] = svd(etawt);
md = min(size(deta));
bigev = find(diag(deta(1:md,1:md))>realsmall);
ueta = ueta(:,bigev);
veta = veta(:,bigev);
deta = deta(bigev,bigev);

if check_exist
    zwt = u2*psi;
    [uz,dz,vz] = svd(zwt);
    md = min(size(dz));
    bigev = find(diag(dz(1:md,1:md))>realsmall);
    uz = uz(:,bigev);
    vz = vz(:,bigev);
    dz = dz(bigev,bigev);

    if isempty(bigev)
        eu(1) = 1;
    else
        disp(norm(uz-ueta*ueta'*uz,'fro'))
        disp(realsmall*n)
        eu(1) = (norm(uz-ueta*ueta'*uz,'fro') < realsmall*n);
    end

    if (~eu(1) & (n_v == -1))
        error('Solution does not exist');
    end
    impact = real(-pi*veta*(deta\ueta')*uz*dz*vz'+psi);
else
    impact = real(-pi*veta*(deta\ueta')*u2*psi+psi);
end

G1 = U*T*spdiags([ones(n-nunstab,1);zeros(nunstab,1)],0,n,n)*U';
G1 = real(G1);

if check_uniq
    [~,deta1,veta1] = svd(u1*pi);
    md = min(size(deta1));
    bigev = find(diag(deta1(1:md,1:md))>realsmall);
    veta1 = veta1(:,bigev);
    if isempty(veta1)
        eu(2) = 1;
    else
        eu(2) = norm(veta1-veta*veta'*veta1,'fro')<realsmall*n;
    end
end
F = u1(:,1:nunstab)'*inv(u1(:,nunstab+1:end)');
C = 1;  	% constant term is not coded yet
