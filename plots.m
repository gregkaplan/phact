% SeHyoun Ahn, Greg Kaplan, Ben Moll, and Tom Winberry
% May 15th, 2016

clear all;
clc;
close all;

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
rrho0 = .01;
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
Ir = 50;									% maximum iterations on steady state interest rate
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
nVars = 2 * I * J - 1 + 5;						% endogenous variables (states + controls)
nEErrors = I * J + 1;						% expectational equations
nShocks = 2;								% number of aggregate shocks

%----------------------------------------------------------------
% Solve for Steady State
%----------------------------------------------------------------

fprintf('Computing steady state...\n')

global IfSS IbSS I0SS varsSS YSS rSS da_tilde A da_stacked

% Compute steady state
t0 = tic;
[rSS,wSS,rKSS,YSS,KSS,mSS,VSS,gSS,cSS,dV_UpwindSS,dVfSS,dVbSS,IfSS,IbSS,I0SS] = computeSteadyState();

tElapsed = toc(t0);

fprintf('...Done!\n')
fprintf('Time to compute steady state: %2.4f seconds\n\n\n',tElapsed)

% Collect endogenous variables in steady state
varsSS = zeros(nVars,1);
varsSS(1:I*J,1) = reshape(VSS,I*J,1);			% value function
ggSS = reshape(gSS,I*J,1);
varsSS(I*J+1:2*I*J-1,1) = ggSS(1:I*J-1,1);		% distribution
varsSS(2*I*J,1) = 0;							% inflation
varsSS(2*I*J+1,1) = YSS;						% output
varsSS(2*I*J+2,1) = mSS;						% marginal costs

% Plot steady state distributions (turn off for speed)
g_a = sum(gSS,2);
g_z = da_tilde'*gSS;

figure
hold on
plot(a / (wSS * (1 - ttau) * z_avg * nSS),gSS(:,cutoff1) / g_z(cutoff1),'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(a / (wSS * (1 - ttau) * z_avg * nSS),gSS(:,cutoff2) / g_z(cutoff2),'linewidth',1.5,'linestyle','-','color',[.8,.24,.03])
xlim([min(a) .2*max(a)])
xlabel('Net worth to average income','interpreter','latex')
ylabel('Mass of households','interpreter','latex')
h = legend('10th percentile income','90th percentile income');
set(h,'interpreter','latex','location','northeast')
set(gcf,'color','w')
grid on
hold off

figure
hold on
f = bar(a / (wSS * (1 - ttau) * z_avg * nSS),sum(gSS,2),'histc')
sh=findall(gcf,'marker','*');delete(sh);
set(f,'FaceColor',[.03,.24,.8],'EdgeColor',[.8,.24,.03]);
text(.5,0.9,['$\leftarrow Pr(a=0)=' num2str((g_a(1).*da_tilde(1))) '$'],'Fontsize',10,'interpreter','latex','Color','k');
xlim([min(a / (wSS * (1 - ttau) * z_avg * nSS)),.2*max(a / (wSS * (1 - ttau) * z_avg * nSS))])
grid on
set(gcf,'color','w')
hold off

popfrac = cumsum(g_a.*da_tilde);
quantile = [0.1 0.01 0.001];
integrand = a.*g_a.*da_tilde;
for p=1:3
    [obj,index(p)] = min(abs(popfrac - (1-quantile(p))));
    wealth_share(p) = sum(integrand(index(p):I))/KSS;
end
fprintf('Fraction with a < 0 = %2.4f\n',100*sum(g_a(1:I_neg).*da_tilde(1:I_neg)))
fprintf('Top 10 percent wealth share = %2.4f\n',100*wealth_share(1))
fprintf('Top 1 percent wealth share = %2.4f\n',100*wealth_share(2))
fprintf('Top .1 percent wealth share = %2.4f\n\n\n',100*wealth_share(3))

% Plot steady state MPCs
tau_MPC = 1;
N_MPC = 21;
dt_MPC = tau_MPC / (N_MPC - 1);

C_stacked = zeros(I*J,N_MPC);
ctilde = zeros(I*J,1);
C_stacked(:,N_MPC) = 0;

B = (1/dt_MPC) * speye(I*J) - A;
c_stacked = reshape(cSS,I*J,1);

for n=N_MPC-1:-1:1
	vec = c_stacked + C_stacked(:,n+1)/dt_MPC;
	C_stacked(:,n) = B\vec;
end

C = reshape(C_stacked(:,1),I,J);
MPC = (C(2:I,:) - C(1:I-1,:)) ./ (aa(2:I,:) - aa(1:I-1,:));

figure
hold on
plot(a(2:I) / (wSS * (1 - ttau) * z_avg * nSS),MPC(:,cutoff1),'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(a(2:I) / (wSS * (1 - ttau) * z_avg * nSS),MPC(:,cutoff2),'linewidth',1.5,'linestyle','-','color',[.8,.24,.03])
xlim([min(a(2:I) / (wSS * (1 - ttau) * z_avg * nSS)) .2*max(a(2:I) / (wSS * (1 - ttau) * z_avg * nSS))])
xlabel('Net worth to average income','interpreter','latex')
ylabel('MPC','interpreter','latex')
h = legend('10th percentile income','90th percentile income');
set(h,'interpreter','latex','location','northeast')
set(gcf,'color','w')
grid on
hold off

%----------------------------------------------------------------
% Plot stationary IRFs
%----------------------------------------------------------------

dt = .005;
T = 100; N_sim = T / dt;
vTime = linspace(1,T,N_sim);
nShock = 2/dt+1;

load irfsSS.mat
load irfsNoMon.mat
load irfsMon.mat

figure
hold on
plot(vTime(nShock:N_sim),100*(vAggConsumptionSS(nShock:N_sim,1) - vAggConsumptionSS(N_sim,1))/vAggConsumptionSS(N_sim,1),'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
xlim([vTime(nShock) vTime(N_sim)])
xlabel('Quarters since shock','interpreter','latex')
ylabel('Percentage deviation from steady state','interpreter','latex')
set(gcf,'color','w')
grid on
hold off

figure
hold on
plot(vTime(nShock:N_sim),100*(vAggConsumptionSS(nShock:N_sim,1) - vAggConsumptionSS(N_sim,1))/vAggConsumptionSS(N_sim,1),'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(vTime(nShock:N_sim),100*(vAggConsumption_mon(nShock:N_sim,1)-vAggConsumption(nShock:N_sim,1))/vAggConsumptionSS(N_sim,1),'linewidth',1.5,'linestyle','-','color',[.8,.24,.03])
xlim([vTime(nShock) vTime(N_sim)])
xlabel('Quarters since shock','interpreter','latex')
ylabel('Percentage deviation from steady state','interpreter','latex')
set(gcf,'color','w')
grid on
hold off