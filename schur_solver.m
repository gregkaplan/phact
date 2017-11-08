function [G1,C,impact,eu,F]=schur_solver(g0,g1,c,psi,pi,continuous,check_exist,check_uniq,n_v)
% Solver based on Schur decomposition (This is based on Sims's gensys, but adjusted
%    to take advantage of special case that applies to heterogeneous agent models).
%
% by SeHyoun Ahn, Dec 2016
%
% REFERENCES:
%    * Sims, Christopher A. "Solving linear rational expectations models."
%         Computational economics 20.1 (2002): 1-20.
%    * Ahn, SeHyoun, Greg Kaplan, Benjamin Moll, Thomas Winberry, and
%         Christian Wolf. "When Inequality Matters for Macro and Macro Matters
%         for Inequality."
%
% PARAMETERS:
%    Check Sims's paper for notation or Our Paper for notation
%
%    XXXXX DO NOT USE DISCRETE TIME VERSION. IT PROBABLY IS WRONG XXXXXX
%    The issue is noted in < https://github.com/gregkaplan/phact/issues/1 >
%
%    continuous = 1 for a continuous time problem (default)
%                0 for a discrete time problem
%
%    check_exist = 1 check for existence of solution (default)
%                 0 existence has been checked before
%    check_uniq = 0 do not check for uniqueness (default)
%                1 check for uniqueness
%    n_v = number of unstable roots
%
% OUPUTS:
%    eu(1) = 1 solution exists
%            0 no stable solution
%           -5 Flag for checking existence is off
%    eu(2) = 1 solution is unique
%            0 solution is not unique
%           -5 Flag for checking uniqueness is off
%
% SYNTAX:
% [G1,C,impact,eu,F]=schur_solver(g0,g1,c,psi,pi,continuous,check_exist,check_uniq,varargin)

%% Set Default Values %%
switch nargin
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
    g_eigs = real(ordeig(T));
    nunstab = sum(g_eigs>0);
    if n_v > -1
        [aux,ind] = sort(g_eigs,'descend');
        locs = ones(n,1);
        locs(ind(1:n_v)) = 0;
        [U,T] = ordschur(U,T,locs);
        if nunstab > n_v
            warning('<schur_solver>: There are more than n_v number of positive eigenvalues with smallest values:');
            disp(aux(n_v+1:nunstab));
        elseif nunstab < n_v
            warning('<schur_solver>: There are less than n_v number of positive eigenvalues:');
            disp(aux(nunstab+1:n_v));
        end
        nunstab = n_v;
    else
        [U,T] = ordschur(U,T,'lhp');
    end
else
    g_eigs = abs(ordeig(T));
    nunstab = sum(g_eigs>1);
    if n_v > -1
        [aux,ind] = sort(g_eigs,'descend');
        locs = ones(n,1);
        locs(ind(1:n_v)) = 0;
        [U,T] = ordschur(U,T,locs);
        if nunstab > n_v
            warning('<schur_solver>: There are more than n_v number of positive eigenvalues with smallest values:');
            disp(aux(n_v+1:nunstab));
        elseif nunstab < n_v
            warning('<schur_solver>: There are less than n_v number of positive eigenvalues:');
            disp(aux(nunstab+1:n_v));
        end
        nunstab = n_v;
    else
	[U,T] = ordschru(U,T,'udi');
    end
end

u1 = U(:,1:n-nunstab)';
u2 = U(:,n-nunstab+1:n)';

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
        eu(1) = (norm(uz-ueta*ueta'*uz,'fro') < realsmall*n);
    end

    if (~eu(1) && (n_v == -1))
        warning('<schur_solver>: Solution does not exist');
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
impact = [F*psi(nunstab+1:end,:);psi(nunstab+1:end,:)];
C = c;  	% constant term is not coded yet
