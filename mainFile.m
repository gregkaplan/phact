%% FIRST-ORDER PERTURBATION SOLUTION OF 2-ASSET ABH ECONOMY

% 2 income states, schurSolver for all, no chance of death

%% HOUSEKEEPING

clear all
close all
clc

currentFolder = pwd;
addpath(strcat(currentFolder,'/Auxiliary_Functions'))
addpath(strcat(currentFolder,'/Subroutines'))
addpath(strcat(currentFolder,'/Factorization'))

id = 'MATLAB:nearlySingularMatrix';
warning('off',id)

%% SOLUTION APPROACH

reduceDistribution = 1;
reduceV = 1;

%% PARAMETERS

%----------------------------------------------------------------
% Set economic parameters 
%----------------------------------------------------------------

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T ...
    nVars nEErrors raSS wSS KSS uSS cSS dSS VSS gSS dab varsSS

% preferences

gamma       = 1;                 % risk aversion
deathrate	= 0;          % probability of death
rho         = 0.014526858925819 + (1/(4*45) - deathrate); % discount rate

% adjustment cost

chi0 = 0.070046;
chi1 = 0.535791;
chi2 = 0.971050;
a_lb = 0.2407;

% production

alpha     = 1/3;    % capital share
delta     = 0.025;   % depreciation rate
rho_tfp   = 0.95;
sigma_tfp = 0.007;
mean_tfp  = 0;

% monetary policy

r_b = 0.005;                 % saving rate, liquid
borrwedge = 0.019875;
r_b_borr = r_b + borrwedge;  % borrowing rate

% fiscal policy

trans = 0.430147477922432;       % transfer to households 
tau_I = 0.25;                % tax rate on wage income

% deposit share

xi = 0;

% perfect annuity markets

pam = 1;

%----------------------------------------------------------------
% Set grids
%----------------------------------------------------------------

% income states

load income_grid
load income_transition

y = [0,1];

lambda1 = 1/2;						 % expected duration of unemployment is 2 quarters
lambda2 = (lambda1 / (y(2) * .93 - y(1))) ...
	* (y(2) - y(2) * .93);				 % unemployment rate 7%
lambda = [-lambda1, lambda1; lambda2, - lambda2];

N = size(y,2);                          % number of income states
% y = exp(y);

y_dist = stat_dist(lambda');            % stationary distribution

ymean = y * y_dist;
y = y./ymean;                           % ensure mean 1

% liquid assets

a_ = 0;
b_ = -2;

I = 40; % number of grid points
I_neg = I/10; % number of grid points on the negative part of the real line
bmin = b_; % lower bound (usually negative, roughly minus the wage)
bmax = 80; % upper bound

b_pos = linspace(0,1,I-I_neg)'; % set the positive part of the grid
coeff_power = 0.9; power = 8;
b_pos = bmax*((1 - coeff_power) * b_pos + coeff_power * (b_pos.^power));

b_neg = linspace(0,1,I_neg+1)'; % set the negative part of the grid
coeff_power = 0; coeff_lin = 0; power = 1;
b_neg = -bmin*((1 - coeff_power) * b_neg + coeff_power * (b_neg.^power)) + bmin;
b_neg(end) = [];

b = [b_neg;b_pos]; % combine the grid

% illiquid assets

J = 30; % number of grid points
amin = a_; % lower bound (usually 0)
amax = 200; % upper bound

a = linspace(0,1,J)'; % set the grid
coeff_power = 0.9; power = 8;
a = amax*((1 - coeff_power) * a + coeff_power * (a.^power));

%----------------------------------------------------------------
% Set initial guesses 
%----------------------------------------------------------------

KL_0 = 20.1215;

%----------------------------------------------------------------
% Time discretization
%----------------------------------------------------------------

T       = 100;							 % number of quarters in simulation
dt      = 0.1;							 % size of time steps
periods = T / dt;

%----------------------------------------------------------------
% Total number of variables
%----------------------------------------------------------------

nVars    = I*J*N*2 - 1 + 4; % I*J*N entries for V and g, one redundancy in g, and K, r, w, tfp
nEErrors = I*J*N; % expectational errors in the HJB evolution

%% STEP 1: STEADY STATE COMPUTATION

% computation

fprintf('Computing steady state...\n')

t0 = tic;
[raSS,wSS,KSS,uSS,cSS,dSS,VSS,gSS,dab] = computeSteadyState();
tElapsed = toc(t0);

fprintf('Time to compute steady state: %.3g seconds\n\n\n',tElapsed)

% store results

varsSS                      = zeros(I*J*N*2 - 1 + 4,1);
varsSS(1:I*J*N,1)           = reshape(VSS,I*J*N,1);
ggSS                        = reshape(gSS,I*J*N,1);
varsSS(I*J*N+1:2*I*J*N-1,1) = ggSS(1:I*J*N-1);
varsSS(2*I*J*N,1)           = KSS;
varsSS(2*I*J*N+1,1)         = raSS;
varsSS(2*I*J*N+2,1)         = wSS;
varsSS(2*I*J*N+3,1)         = 0; % tfp

% plot asset holdings

plot_dist

%% STEP 2: LINEARIZE MODEL EQUATIONS

% Compute derivatives w.r.t. endogenous variables
vars = zeros(nVars + nVars + nEErrors + 1,1); % values differentiating at: V, g, K, r, w, Z, Vdot, gdot, empty 3 x 1, Zdot, eta, W
vars = myAD(vars);		% converts to dual numbers, required for differentiation
derivativesIntermediate = equilibriumConditions(vars);	% evaluates function and derivatives
derivs = getderivs(derivativesIntermediate);	% extracts only derivatives

% Unpackage derivatives
mVarsDerivs = derivs(:,1:nVars);
mVarsDotDerivs = derivs(:,nVars+1:2*nVars);
mEErrorsDerivs = derivs(:,2*nVars+1:2*nVars+nEErrors);
mShocksDerivs = derivs(:,2*nVars+nEErrors+1);

%% STEP 3: SOLVE LINEAR SYSTEM

% define matrices

g0 = mVarsDotDerivs;
g1 = -mVarsDerivs;
c = sparse(nVars,1);
psi = -mShocksDerivs;
pi = -mEErrorsDerivs;

% some useful quantities

n_v = I*J*N;
n_g = I*J*N-1;
n_p = 3;
n_Z = 1;

if reduceDistribution == 0

    % invert matrix

    [base,inv_base,g1,c,pi,psi] = cleanG0Sparse(g0,g1,c,pi,psi);
    
else
    
    base = speye(nVars);
    inv_base = speye(nVars);
    
end

% distribution reduction

if reduceDistribution == 1

    [state_red,inv_state_red,n_g] = stateSpaceReduction(g0,g1,n_v,n_g,n_p,n_Z);
    [g1,psi,pi,c] = changeBasis(state_red,inv_state_red,g1,psi,pi,c);
    
elseif reduceDistribution == 0

    state_red = speye(nVars);
    inv_state_red = speye(nVars);
    
end

% VF reduction

if reduceV == 1
    
    c_a_power = 5;
    n_a_knots = 5; % number of knots

    c_b_power = 5; 
    n_b_knots = 5; % number of knots
    
    b_knots = linspace(0,1,n_b_knots)';
    b_knots = (exp(c_b_power*b_knots)-1)/(exp(c_b_power)-1);
    
    b_knots = min(b) + (max(b) - min(b)) * b_knots;

    a_knots = linspace(0,1,n_a_knots)';
    a_knots = (exp(c_a_power*a_knots)-1)/(exp(c_a_power)-1);
    
    a_knots = min(a) + (max(a) - min(a)) * a_knots;
    
    [from_spline_small,to_spline_small] = quad2DsplineND(b,a,b_knots,a_knots);

    n_splined = size(from_spline_small,2)*N;
    n_rest = n_g + n_Z;

    to_spline_big = kron(speye(N),to_spline_small);
    to_spline = spdiags(ones(n_splined+n_rest,1),n_v-n_splined,n_rest+n_splined,n_v+n_rest);
    to_spline(1:n_splined,1:n_v) = to_spline_big;

    from_spline_big = kron(speye(N),from_spline_small);
    from_spline = spdiags(ones(n_splined+n_rest,1),n_splined-n_v,n_v+n_rest,n_rest+n_splined);
    from_spline(1:n_v,1:n_splined) = from_spline_big;

    [g1,psi,pi,c]=changeBasis(to_spline,from_spline,g1,psi,pi,c);
    
elseif reduceV == 0
    
    if reduceDistribution == 1

        from_spline = speye(n_g + n_v + n_Z);
        to_spline = speye(n_g + n_v + n_Z);
        n_splined = n_v;
    
    elseif reduceDistribution == 0
        
        from_spline = speye(n_g + n_v + n_Z + n_p);
        to_spline = speye(n_g + n_v + n_Z + n_p);
        n_splined = n_v;
        
    end
    
end

% solve system

[G1,~,impact,eu] = schurSolver(speye(n_splined+n_g+n_Z),g1,c,psi,pi,1,1,1);
n_rest=n_g+n_Z;
disp(eu)

tSolver = toc(t0);

%% STEP 4: SIMULATE SYSTEM

% response to aggregate technology shock

N = T / dt;
n_v = n_splined;
vTime = linspace(1,T,N);
v_small = zeros(n_rest+n_v,N+1);
vAggregateShock = zeros(N,1);
vAggregateShock(1:1/dt,1) = 1;

t0=tic;

GG1 = (speye(size(G1,1)) - G1 * dt)^(-1);

for n = 1 : N
	v_small(:,n+1) = GG1 * (v_small(:,n) + (dt^(1/2)) * impact * vAggregateShock(n,1)');
end

if reduceDistribution == 0 && reduceV == 1
    vVarsSeries = base*from_spline*v_small;
elseif reduceDistribution == 0 && reduceV == 0
    vVarsSeries = base*v_small;
else
    vVarsSeries = inv_state_red*from_spline*v_small;
end

fprintf('Time to simulate Reduced Model: %.3g seconds\n',toc(t0));

%% PLOT IRFs

figure(2)
hold on
plot(vTime(1/dt:end),100*vVarsSeries(4*I*J+3,1/dt+1:end),'k-','linewidth',1.5)
plot(vTime(1/dt:end),100*(log(vVarsSeries(4*I*J,1/dt+1:end) + KSS) - log(KSS)),'b-','linewidth',1.5)
legend('TFP','Capital')
title('Impulse response to positive TFP shock')
hold off

%% SAVE RESULTS

time_IRF = vTime(1/dt:end);
K_IRF = 100*(log(vVarsSeries(4*I*J,1/dt+1:end) + KSS) - log(KSS));

if reduceDistribution == 0 && reduceV == 0
    save IRF_v0g0_v3 time_IRF K_IRF
elseif reduceDistribution == 1 && reduceV == 0
    save IRF_v0g1_v3 time_IRF K_IRF
elseif reduceDistribution == 0 && reduceV == 1
    save IRF_v1g0_v3 time_IRF K_IRF
elseif reduceDistribution == 1 && reduceV == 1
    save IRF_v1g1_v3 time_IRF K_IRF
end

%% COMPARISON PLOT

clc
clear all

load IRF_v0g0_v3
K_IRF_v0g0 = K_IRF;

load IRF_v0g1_v3
K_IRF_v0g1 = K_IRF;

load IRF_v1g0_v3
K_IRF_v1g0 = K_IRF;

load IRF_v1g1_v3
K_IRF_v1g1 = K_IRF;

figure(3)
hold on
plot(time_IRF,K_IRF_v0g0,'b-','linewidth',1.5)
plot(time_IRF,K_IRF_v1g0,'g-','linewidth',1.5)
plot(time_IRF,K_IRF_v0g1,'y-','linewidth',1.5)
plot(time_IRF,K_IRF_v1g1,'r-','linewidth',1.5)
legend('v0g0','v1g0','v0g1','v1g1')
title('Impulse response to positive TFP shock')
hold off