function [G1,C,impact,eu]=schurSolver(g0,g1,c,psi,pi,continuous,check_exist,check_uniq)
% Solver based on Schur decomposition
%
% by SeHyoun Ahn, Dec 2016
%
%% PARAMETERS:
% continuous= 1 for a continuous time problem (default)
%             0 for a discrete time problem
% check_exist= 1 check for existence of solution (default)
%              0 existence has been checked before
% check_uniq= 0 do not check for uniqueness (default)
%             1 check for uniqueness
%
%% OUPUTS:
% eu(1)= 1 solution exists
%        0 no stable solution
%       -5 Flag for checking existence is off
% eu(2)= 1 solution is unique
%        0 solution is not unique
%       -5 Flag for checking uniqueness is off

%% Set Default Values %%
switch nargin
    case 8
    case 7
        check_uniq=0;
    case 6
        check_uniq=0;
        check_exist=1;
    case 5
        check_uniq=0;
        check_exist=1;
        continuous=1;
    case 1:4
        error('Insufficient Number of Input Variables');
end
realsmall=sqrt(eps)*10;
eu=[-5;-5];
n=size(g1,1);

%% Schur Decomposition %%
[U,T]=schur(full(g1),'real');
if continuous
    [U,T]=ordschur(U,T,'lhp');
    g_eigs=ordeig(T);
    nunstab=sum(g_eigs>0);
else
    [U,T]=ordschru(U,T,'udi');
    g_eigs=abs(ordeig(T));
    nunstab=sum(g_eigs<1);
end

u2=U(:,n-nunstab+1:n)';
u1=U(:,1:n-nunstab)';

etawt=u2*pi;
[ueta,deta,veta]=svd(etawt);
md=min(size(deta));
bigev=find(diag(deta(1:md,1:md))>realsmall);
ueta=ueta(:,bigev);
veta=veta(:,bigev);
deta=deta(bigev,bigev);

if check_exist
    zwt=u2*psi;
    [uz,dz,vz]=svd(zwt);
    md=min(size(dz));
    bigev=find(diag(dz(1:md,1:md))>realsmall);
    uz=uz(:,bigev);
    vz=vz(:,bigev);
    dz=dz(bigev,bigev);

    if isempty(bigev)
        eu(1)=1;
    else
        eu(1)=(norm(uz-ueta*ueta'*uz,'fro') < realsmall*n);
    end

    if ~eu(1)
        error('Solution does not exist');
    end
    impact=real(-pi*veta*(deta\ueta')*uz*dz*vz'+psi);
else
    impact=real(-pi*veta*(deta\ueta')*u2*psi+psi);
end

G1=U*T*spdiags([ones(n-nunstab,1);zeros(nunstab,1)],0,n,n)*U';
G1=real(G1);
    
if check_uniq
    [~,deta1,veta1]=svd(u1*pi);
    md=min(size(deta1));
    bigev=find(diag(deta1(1:md,1:md))>realsmall);
    veta1=veta1(:,bigev);
    if isempty(veta1)
        eu(2)=1;
    else
        eu(2)=norm(veta1-veta*veta'*veta1,'fro')<realsmall*n;
    end
end
C=1; 	% constant term is not coded yet
