function [G1,C,impact,eu]=phact_solver(g0,g1,c,psi,pi,clean,continuous,check_exist,check_uniq)
%% FLAGS:
% clean= 1 g0 is not an identity matrix (default)
%        0 g0 is already an identity matrix
% continuous= 1 for a continuous time problem (default)
%             0 for a discrete time problem
% check_exist= 1 check for existence of solution (default)
%              0 existence has been checked before
% check_uniq= 0 do not check for uniqueness (default)
%             1 check for uniqueness

%% OUPUT:
% eu(1)= 1 solution exists
%        0 no stable solution
%       -5 Flag for checking existence is off
% eu(2)= 1 solution is unique
%        0 solution is not unique
%       -5 Flag for checking uniqueness is off


%% Set Default Values %%
global nunstab
switch nargin
    case 9
    case 5
        check_uniq=0;
        check_exist=1;
        continuous=1;
        clean=1;
    case 8
        check_uniq=0;
    case 7
        check_uniq=0;
        check_exist=1;
    case 6
        check_uniq=0;
        check_exist=1;
        continuous=1;
    case 1:4
        error('Insufficient Number of Input Variables');
end
realsmall=sqrt(eps)*10;
eu=[-5;-5];

%% Get Input into Schur Form %%
if clean
    %disp('Converting to Reduced Form');
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

    lastwarn('');
    %g1=g0(keep,keep)\g1(keep,:)*base_change;
    %psi=g0(keep,keep)\psi(keep,:);
    %pi=g0(keep,keep)\pi(keep,:);
    %c=g0(keep,keep)\c(keep,:);
    [~, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:singularMatrix')
        error('Wrong Form. Try running Gensys');
    end
else
    lastwarn('');
    g1=g0\g1;
    [~, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:singularMatrix')
        error('g0 is not invertible. Try running with <clean=1>');
    end
end
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
        %error('Solution does not exist');
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

if clean
    %tmp=speye(size(base_change,1));
    G1=base*G1*base';%*tmp(keep,:);
    impact=base*impact;
end
C=1;