% SeHyoun Ahn, Greg Kaplan, Ben Moll, and Tom Winberry
% May 15th, 2016

clear;
clc;
close all;

%% SOLUTION APPROACH

reduceDistribution = 1;
reduceV = 1;
%----------------------------------------------------------------
% Set economic parameters 
%----------------------------------------------------------------

% Preferences
global ggamma rrho pphi ddeath cchi nSS rrho0 rrhoMin rrhoMax
ggamma = 2;									% coefficient of relative risk aversion
pphi = 2;									% inverse frisch elasticity
ddeath = 1 / (4 * 45);						% death rate
nSS = 1 / 3;								% steady state hours
%rrho0 = .06 / 4;							% initial guess of time preference (for calibration)
rrho0 = .02;
rrhoMin = .01 / 4;							% lower bound of time preference (for calibration)
rrhoMax = .15 / 4;							% upper bound of time preference (for calibration)

% Technology
global aalpha ddelta
aalpha = 1 / 3;								% capital coefficient
ddelta = .025;								% depreciation rate

% Financial wedge
global wedge
%wedge = .06/4;								% wedge between borrowing and lending
wedge = 0;

% Tax code
global ttransfer ttau
ttau = 0.25;								% marginal tax rate

% Firm parameters
global eepsilon ttheta
eepsilon = 10;								% elasticity of demand
ttheta = 100;								% price stickiness

% Monetary policy
global pphiInflation pphiOutput
pphiInflation = 1.5;						% coefficient on inflation
pphiOutput = .5;							% coefficient on output gap

% Aggregate shocks
global nnuTFP ssigmaTFP nnuMon ssigmaMon
ssigmaTFP = .007; nnuTFP = 1 - .95;			% total factor productivity
ssigmaMon = .0025; nnuMon = 1 - .05;		% monetary policy shocks

%----------------------------------------------------------------
% Income shock process
%----------------------------------------------------------------

global z J z_avg g_z

% Load grid and transition probabilities from calibration
load ygrid_combined.txt
load ymarkov_combined.txt

% Create grid of income shocks
z = exp(ygrid_combined);
[J,~] = size(z);

% Compute stationary income distribution
AT = ymarkov_combined';
g_z = ones(J,1) / J;
for n = 1 : 50
	% Implicit method for updating distribution
	g_z_new = (speye(J) - AT * 1000) \ g_z;
	diff = max(abs(g_z_new - g_z));
	if diff < 10e-6
		break
	end
	g_z = g_z_new;
end

% Save results into stationary distribution g_y
z_avg = sum(z' * g_z);
z = z';

% Compute 40th percentile
global z_cutoff
zPercentiles = cumsum(g_z);
z_cutoff_index = knnsearch(zPercentiles,.4);
z_cutoff = z(z_cutoff_index);
cutoff1 = knnsearch(zPercentiles,.1);
cutoff2 = knnsearch(zPercentiles,.9);

%----------------------------------------------------------------
% Set other approximation parameters
%----------------------------------------------------------------

% Wealth grid (will build the grid in steady state routine)
global I I_neg amin amax a aa
I = 40; I_neg = 0;						% size of asset grid space

% Expand labor productivity grid
global zz zzz
zz = ones(I,1) * z;
zzz = reshape(zz,I * J,1);

% Matrix which will be helpful for later computations
global Aswitch
Aswitch = kron(ymarkov_combined,speye(I));

% Steady state computations and calibration
global rmin rmax r0 maxit crit Delta Ir crit_S Irrho crit_objective
rmin = .0001;								% lower bound for steady state interest rate
rmax = rrho0;								% upper bound for steady state interest rate
r0 = .005;									% initial guess for steady state interest rate
maxit = 100;								% maximum iterations on steady state HJB
crit = 1e-6;								% error criterion for steady state value function convergence
Delta = 1000;								% update size for implicit scheme on steady state HJB
Ir = 500;									% maximum iterations on steady state interest rate
crit_S = 1e-5;								% error criterion for steady state interest rate
Irrho = 1;									% maximum iterations on steady state calibbration
crit_objective = 1e-3;						% error criterion for steady state calibration

% Time discretization for simulation
global T dt N_sim 
T = 40;										% number of quarters for impulse responses
dt = .1;									% size of time step
N_sim = T / dt;								% number of time periods to simulate

% Number of variables in the system
global nVars nEErrors nShocks
nShocks = 1;								% number of aggregate shocks
nVars = 2 * I * J - 1 + 3 + nShocks;			% endogenous variables (states + controls)
nEErrors = I * J;						% expectational equations

%----------------------------------------------------------------
% Solve for Steady State
%----------------------------------------------------------------
fprintf('Computing steady state...\n')

global IfSS IbSS I0SS varsSS YSS rSS da_tilde A da_stacked

% Compute steady state
t0 = tic;
[rSS,wSS,YSS,KSS,VSS,gSS,cSS,dV_UpwindSS,dVfSS,dVbSS,IfSS,IbSS,I0SS] = computeSteadyState_flex();

tElapsed = toc(t0);

fprintf('...Done!\n')
fprintf('Time to compute steady state: %2.4f seconds\n\n\n',tElapsed)


% Collect endogenous variables in steady state
varsSS = zeros(nVars,1);
varsSS(1:I*J,1) = reshape(VSS,I*J,1);			% value function
ggSS = reshape(gSS,I*J,1);
varsSS(I*J+1:2*I*J-1,1) = ggSS(1:I*J-1,1);		% distribution
varsSS(2*I*J,1)           = KSS;
varsSS(2*I*J+1,1)         = rSS;
varsSS(2*I*J+2,1)         = wSS;
varsSS(2*I*J+3,1)         = 0; % tfp

%----------------------------------------------------------------
% Compute derivatives of equilibrium conditions
%----------------------------------------------------------------

fprintf('Taking derivatives of equilibrium conditions...\n')

t0 = tic;

% Compute derivatives
vars = zeros(2*nVars+nEErrors+nShocks,1);		% derivaties w.r.t current and lagged variables + expectational errors + shocks
vars = myAD(vars);								% convert to dual numbers, for automatic differentiation package
derivativesIntermediate = equilibriumConditions_flex_big(vars);	% evaluates function and derivatives
derivs = getderivs(derivativesIntermediate);	% extracts derivatives

% Unpackage derivatives
mVarsDerivs = derivs(:,1:nVars);
mVarsDotDerivs = derivs(:,nVars+1:2*nVars);
mEErrorsDerivs = derivs(:,2*nVars+1:2*nVars+nEErrors);
mShocksDerivs = derivs(:,2*nVars+nEErrors+1);

fprintf('...Done!\n')
fprintf('Time to compute derivatives: %2.4f seconds\n\n\n',toc(t0))

% define matrices

g0 = mVarsDotDerivs;
g1 = -mVarsDerivs;
constant = sparse(nVars,1);
psi = -mShocksDerivs;
pi = -mEErrorsDerivs;

% some useful quantities

n_v = I*J;
n_g = I*J-1;
n_p = 3;
n_Z = 1;

% % clean G0
% 
% [base,inv_base,g1,constant,pi,psi] = cleanG0Sparse(g0,g1,constant,pi,psi);
% [nVarsCleaned,~] = size(g1);
% g0 = speye(nVarsCleaned);

% distribution reduction

if reduceDistribution == 1

    [state_red,inv_state_red,n_g] = stateSpaceReduction(g0,g1,n_v,n_g,n_p,n_Z);
    [g1,psi,pi,constant] = changeBasis(state_red,inv_state_red,g1,psi,pi,constant);
    
elseif reduceDistribution == 0

    state_red = speye(nVars);
    inv_state_red = speye(nVars);
    
end

% value function reduction

if reduceV == 1
    
    n_knots = 12;
    c_power = 7;
    x = a';
    n_post = 33;	        % This is from two income states
    n_prior = 1;
    % Create knot points for spline (the knot points are not uniformly spaced)
    knots = linspace(amin,amax,n_knots-1)';
    knots = (amax-amin)/(2^c_power-1)*((knots-amin)/(amax-amin)+1).^c_power+amin-(amax-amin)/(2^c_power-1);
    [from_spline_small, to_spline_small] = createQuadSplineNDa(x,knots,n_prior,n_post,n_g);
    
    n_splined = size(from_spline_small,2)*J;
    if reduceDistribution == 1
        n_rest = n_g + n_Z;
    elseif reduceDistribution == 0
        n_rest = n_g + n_Z + n_p;
    end

    to_spline_big = kron(speye(J),to_spline_small);
    to_spline = spdiags(ones(n_splined+n_rest,1),n_v-n_splined,n_rest+n_splined,n_v+n_rest);
    to_spline(1:n_splined,1:n_v) = to_spline_big;

    from_spline_big = kron(speye(J),from_spline_small);
    from_spline = spdiags(ones(n_splined+n_rest,1),n_splined-n_v,n_v+n_rest,n_rest+n_splined);
    from_spline(1:n_v,1:n_splined) = from_spline_big;
    
    if reduceDistribution == 1
        
        [g1,psi,pi,constant]=changeBasis(to_spline,from_spline,g1,psi,pi,constant);
        
    elseif reduceDistribution == 0
    
        [g1,psi,pi,constant,g0]=changeBasis(to_spline,from_spline,g1,psi,pi,constant,g0);
        
    end
    
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

%----------------------------------------------------------------
% Solve problem
%----------------------------------------------------------------

t0 = tic;

fprintf('Solving reduced linear system...\n')

if reduceDistribution == 0
    [G1,C,impact,eu] = phact_solver(g0,g1,constant,psi,pi,1,1,1,1);
    n_rest = n_g + n_Z + n_p;
else
    [G1,C,impact,eu]=phact_solver(speye(n_splined+n_g+n_Z),g1,constant,psi,pi,1,1,1,1);
    n_rest=n_g+n_Z;
    disp(eu)
end

fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
fprintf('Time to solve reduced system: %2.4f seconds\n\n\n',toc(t0))

%----------------------------------------------------------------
% Compute impulse responses to productivity shock
%----------------------------------------------------------------

dt = .01;

T = 100; N_sim = T / dt;

t0 = tic;
vTime = linspace(1,T,N_sim);

fprintf('Simulating response to TFP shock...\n')

% Initialize vector of endogenous variables
v_small = zeros(n_rest+n_splined,N_sim+1);
vAggregateShock = zeros(N_sim,1);
vAggregateShock(1:1/dt,1) = 1;
nShock = 1/dt+1;

% Preallocate other vectors to store series of interest
vAggConsumption = zeros(N_sim,1);
vAggConsumptionPoor = zeros(N_sim,1);
vAggConsumptionRich = zeros(N_sim,1);
vAggOutput = zeros(N_sim,1);
vAggInvestment = zeros(N_sim,1);
vR_real = zeros(N_sim,1);
vTFP = zeros(N_sim,1);
vK = zeros(N_sim,1);

GG1 = (speye(size(G1,1)) - G1 * dt)^(-1);

% Simulate
for n = 1 : N_sim
    
    if mod(n,100) == 0
        disp(n)
    end

	% Compute evolution of endogeneous variables
    v_small(:,n+1) = GG1 * (v_small(:,n) + (dt^(1/2)) * impact * vAggregateShock(n,1)');
end

vVarsSim = inv_state_red*from_spline*v_small;

for n = 1 : N_sim-1

	%----------------------------------------------------------------
	% Compute decisions given known objects
	%----------------------------------------------------------------

    % Extract objects
	vars = vVarsSim(:,n+1);
	V = reshape(vars(1:I*J,1) + varsSS(1:I*J,1),I,J);
	gg = vars(I*J+1:2*I*J-1) + varsSS(I*J+1:2*I*J-1);
	g_end = (1 / da_stacked(I*J,1)) * (1 - sum(gg .* da_stacked(1:I*J-1)));
	g = [gg;g_end];
	logAggregateTFP = vars(2*I*J+3,1);
	
	%%%
	% Compute auxiliary variables from known relationships
	%%%
	
	K = (reshape(aa,I*J,1) .* g)'*da_stacked;
	N = (((1 - ttau) / cchi) * (z_avg ^ pphi) * (1 - aalpha) * exp(logAggregateTFP) * (K ^ aalpha)) ^ (1 / (pphi + aalpha));
	output = exp(logAggregateTFP) * (K ^ aalpha) * (N ^ (1 - aalpha));
	r_real = aalpha * output / K - ddelta;
	w = (1 - aalpha) * output / N;
	nHours = N / z_avg;

	% Adjust rates of return by death
	rrho_effective = rrho + ddeath;
	r = ddeath + r_real .* (aa >= 0) + (r_real + wedge) .* (aa < 0);
	
	%%% 
	% Solve HJB to get consumption and savings
	%%%
	
	c0 = w * (1 - ttau) * nHours * zz + r .* aa + ttransfer;
	disutil = cchi * zz .* ((nHours ^ (1 + pphi)) / (1 + pphi));

	% Compute forward difference
	dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) ./ (aa(2:I,:) - aa(1:I-1,:));
	dVf(I,:) = (c0(I,:) - disutil(I,:)) .^ (-ggamma);

	% Compute backward difference
	dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) ./ (aa(2:I,:) - aa(1:I-1,:));
	dVb(1,:) = (c0(1,:) - disutil(1,:)) .^ (-ggamma);

	% Consumption and savings with forward difference
	cf = dVf .^ (-1 / ggamma) + disutil;
	sf = c0 - cf;

	% Consumption and savings with backward difference
	cb = dVb .^ (-1 / ggamma) + disutil;
	sb = c0 - cb;

	% Consumption with no drift
	dV0 = (c0 - disutil) .^ (-ggamma);

	% Consumption choice with upwind differences
	c = IfSS .* cf + IbSS .* cb + I0SS .* c0;
	u = ((c - disutil) .^ (1 - ggamma)) ./ (1 - ggamma);
	savings = c0 - c;
	c_poor = c .* (zz <= cutoff1);
	c_rich = c .* (zz >= cutoff2);

	%%%
	% Compute aggregate variables of interest
	%%%
	
	vAggConsumption(n,1) = (reshape(c,I*J,1) .* g)'*da_stacked;
	vAggConsumptionPoor(n,1) = (reshape(c_poor,I*J,1) .* g)'*da_stacked;
	vAggConsumptionRich(n,1) = (reshape(c_rich,I*J,1) .* g)'*da_stacked;
	vAggInvestment(n,1) = (reshape(savings,I*J,1) .* g)'*da_stacked;
	vAggOutput(n,1) = output;
	vR_real(n,1) = r_real;
	vTFP(n,1) = logAggregateTFP;
    vK(n,1) = K;
	%}
end
	
tSimulate = toc(t0);
fprintf('Time to simulate model: %.3g seconds\n\n',tSimulate);

figure(1);
hold on;
plot(vTime(1:end-1),vK(1:end-1)-KSS,'--')

% %% SAVE RESULTS
% 
% if reduceDistribution == 0 && reduceV == 0
%     
%     vKK_v0g0 = vVarsSim(2640,2:end);
%     save vKK_v0g0 vKK_v0g0
%     
% elseif reduceDistribution == 1 && reduceV == 0
%     
%     vKK_v0g1 = vVarsSim(2640,2:end);
%     save vKK_v0g1 vKK_v0g1
%     
% elseif reduceDistribution == 0 && reduceV == 1
%     
%     vKK_v1g0 = vVarsSim(2640,2:end);
%     save vKK_v1g0 vKK_v1g0
%     
% elseif reduceDistribution == 1 && reduceV == 1
%     
%     vKK_v1g1 = vVarsSim(2640,2:end);
%     save vKK_v1g1 vKK_v1g1
%     
% end
% 
% save vTime vTime
% 
% %% PLOT RESULTS
% 
% clc
% clear all
% 
% load vTime
% load vKK_v0g0
% load vKK_v1g0
% load vKK_v0g1
% load vKK_v1g1
% 
% figure(1)
% hold on
% plot(vTime(1:end-1),vKK_v0g0(1:end-1),'b-','linewidth',1.5)
% plot(vTime(1:end-1),vKK_v1g0(1:end-1),'g-','linewidth',1.5)
% plot(vTime(1:end-1),vKK_v0g1(1:end-1),'y-','linewidth',1.5)
% plot(vTime(1:end-1),vKK_v1g1(1:end-1),'r-','linewidth',1.5)
% legend('v0g0','v1g0','v0g1','v1g1')
% title('Impulse response to positive TFP shock')
% hold off