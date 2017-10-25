%% One-Asset HANK model
% This program solves one-asset HANK model using the pertrubation method outlined
% in Ahn et al. (2017). The model setup is explained in more detail in the
% documentation file provided at \</examples/one_asset_HANK/docs/one_asset_HANK.pdf\>.
%
% The solution will also follow the standard steps
%
% # Solve the steady-state problem
% # Linearize Model Equations
% # Solve out Static Constraints or Reduce Dimensionality
% # Solve Linear System
% # Compute Impulse Response Functions
% # Internal Consistency Check
%
% Estimated Runtime: 1 ~ 5 seconds
%
% REFERENCES:
%
% * Ahn, SeHyoun, et al. "When Inequality Matters for Macro and Macro Matters for
%   Inequality." NBER Macroeconomics Annual 2017, volume 32. University of
%   Chicago Press, 2017.
% * Kaplan, Greg, Benjamin Moll, and Giovanni L. Violante. Monetary policy
%   according to HANK. No. w21897. National Bureau of Economic Research, 2016.
%
% REQUIRES:
%
% * auto diff toolbox: <https://github.com/sehyoun/MATLABAutoDiff>
% * phact toolbox: <https://github.com/gregkaplan/phact>
% * <compute_steady_state.m>
% * <equilibrium_conditions.m>
% * <set_parameters.m>
% * \< plot_IRFs.m \>
%


%%
% Setup the toolbox
% Need to include folders containing the files
%addpath('/path/to/PHACT');
%addpath('/path/to/AutoDiff');

addpath('/home/sehyoun/Dropbox/1Packages/PHACT');
addpath('/home/sehyoun/Dropbox/1Packages/AutoDiff');

% Turn off the warning message from auto diff
% You should read the warning once, but turn it off after reading it.
warning('off','AutoDiff:maxmin')

%%
% Set options for this example run
ReduceDistribution = 1;  % 1 for state space reduction 0 for not
reduceV = 1;             % 1 for value function reduction 0 for not
ReduceDist_hor = 20;     % Dimensionality of the Krylov subspace
DisplayLev = 1;          % Determines verbosity of steady state calculation
check_consistency = 1;   % Runs Step 6: Internal consistency check

%% Step 0: Set parameters
% The script sets up parameters relevant for the model
%    Economic parameters, approximation parameters, and grid parameters are
%    defined in the script.
set_parameters;

n_v = I*J + 1;          % number of jump variables (value function + inflation)
n_g = I*J;              % number of endogeneous state variables (distribution + monetary policy)
n_p = 5;                % number of static relations: bond-market clearing, labor market clearing, consumption, output, total assets
n_shocks = 1;           % only monetary policy shock is considered
nEErrors = n_v;
nVars = n_v + n_g + n_p;

%% Step 1: Solve for the Steady State
% Any methods can be used to solved for the steady state. In particular, example
%    codes can be found at \<http://www.princeton.edu/~moll/HACTproject.html\>.
fprintf('Computing steady state...\n');
t0 = tic;
compute_steady_state;
fprintf('Time to compute steady state: %2.4f seconds\n\n\n',toc(t0));


%% Step 2: Linearize Model Equations
% For computing derivatives, the codes written for solving for the
%    steady-state can be used almost verbatim using automatic
%    differentiation toolbox as long as only the functions supported by
%    automatic differentation are used. For list of supported functions and
%    documentation of relevant syntax check
%    <https://github.com/sehyoun/MATLABAutoDiff>. Example usage/syntax of
%    automatic differentiation can be found at
%    <https://sehyoun.com/EXAMPLE_AutoDiff_Syntax.html>
fprintf('Taking derivatives of equilibrium conditions...\n');
t0 = tic;

% Prepare automatic differentiation
vars = zeros(nVars + nVars + nEErrors + n_shocks,1);
vars = myAD(vars);

% Evaluate derivatives
equilibrium_conditions;

% Extract derivative values
derivs = getderivs(v_residual);

t_derivs = toc(t0);
fprintf('Time to compute derivatives: %2.4f seconds\n\n\n',t_derivs);
if t_derivs > 1
    warning('If you do not compile mex/C files for the automatic differentiation, matrix vector multiplication will be slow');
    disp('Press any key to continue...');
    pause();
end

%% Step 3: Solve out Static Constraints or Reduce the Model
% Extract derivatives
g1 = -derivs(:,1:nVars);
g0 = derivs(:,nVars+1:2*nVars);
pi = -derivs(:,2*nVars+1:2*nVars+nEErrors);
psi = -derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+n_shocks);
constant = sparse(nVars,1);

% State Variables
if ReduceDistribution	== 1
    % Reduce model
	fprintf('Reducing distribution ...\n');
    [state_red,inv_state_red,n_g_red] = krylov_reduction(g0,g1,n_v,n_g,ReduceDist_hor);
    [g1,psi,pi,constant,g0] = change_basis(state_red,inv_state_red,g1,psi,pi,constant,g0);
else
    % Solve out static constraints
    fprintf('Solving Out Static Constraints ...\n');
    [state_red,inv_state_red,g0,g1,constant,pi,psi] = clean_G0_sparse(g0,g1,constant,pi,psi);
    n_g_red = n_g;
end

% Jump Variables
if reduceV == 1
    % Reduce dimensionality of value function using splines
    n_knots = 12;
    c_power = 1;
    x = a';
    n_post = size(z,2);
    n_prior = 1;

    % Create knot points for spline (the knot points are not uniformly spaced)
    knots = linspace(amin,amax,n_knots-1)';
    knots = (amax-amin)/(2^c_power-1)*((knots-amin)/(amax-amin)+1).^c_power+ amin-(amax-amin)/(2^c_power-1);

    % Function calls to create basis reduction
    [from_spline, to_spline] = oneDquad_spline(x,knots);
    [from_spline, to_spline] = extend_to_nd(from_spline,to_spline,n_prior,n_post);
    from_spline(end+1,end+1) = 1;
    to_spline(end+1,end+1) = 1;
    n_splined = size(from_spline,2);
    [from_spline, to_spline] = projection_for_subset(from_spline,to_spline,0,n_g_red);

    % Reduce the decision vector
    [g1,psi,~,constant,g0] = change_basis(to_spline,from_spline,g1,psi,pi,constant,g0);
    pi = to_spline * pi * from_spline(1:n_v,1:n_splined);

else
    % Create identity matrix for code reuse below
    from_spline = speye(n_g_red + n_v);
    to_spline = speye(n_g_red + n_v);
    n_splined = n_v;
end


%% Step 4: Solve Linear Systems
t0 = tic;
fprintf('Solving linear system...\n');

% Note that I (SeHyoun) will probably swap out some parts of schur_solver in the
%    near future, so it might have a different interface. Since the codebase is
%    new, there will be interative updates to simplify interface and syntax.
%    Underlying math will stay the same, but interfact may change with updates,
%    so one should note the release number of the codebase one is using.
[G1, ~, impact, eu, F] = schur_solver(g0,g1,c,psi,pi,1,1,1,n_splined);
fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0));


%% Step 5: Simulate Impulse Response Functions
fprintf('Simulating Model...\n');
t0 = tic;

T = 100;
N = 400;
dt = T/N;
vAggregateShock	= zeros(1,N);
vAggregateShock(1) = 1/sqrt(dt);
trans_mat = inv_state_red*from_spline;
[simulated,vTime] = simulate(G1,impact,T,N,vAggregateShock,'implicit',trans_mat,[n_v,n_v+n_g:n_v+n_g+5]);

fprintf('...Done!\n');
fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0));

inflation = simulated(1,:)';
monetary_shock = simulated(2,:)';
consumption = (simulated(4,:)')/vars_SS(n_v+n_g+3);
Y = (simulated(6,:)')/vars_SS(n_v+n_g+4);
lab_sup = (simulated(4,:)')/vars_SS(n_v+n_g+2);
wage = simulated(3,:)'/vars_SS(n_v+n_g+1);


%% Step 6: Internal Consistency Check
% A different internal consistency check will be implemented and updated in the
%    future. This function should only be taken as a sanity check.

if check_consistency
    g1 = -derivs(:,1:nVars);
    psi = -derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+n_shocks);
    from_red = inv_state_red * from_spline;
    to_red = to_spline * state_red;
    [epsilon] = internal_consistency_check(G1,impact,n_g_red,from_red,to_red,g1,psi,F,n_v,n_g,1000,vars_SS,1,0.07);
end

%% (optional) Step 7: Plot relevant IRFs
line_style = '-';
color = 'blue';
line_style = '--';
color = 'red';

figure(1);
subplot(2,3,1);
hold on;
plot(vTime,monetary_shock,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Monetary Policy Shock','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(2,3,2);
hold on;
plot(vTime,inflation,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Inflation','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(2,3,3);
hold on;
plot(vTime,consumption,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Consumption','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(2,3,4);
hold on;
plot(vTime,Y,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('GDP','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(2,3,5);
hold on;
plot(vTime,lab_sup,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Labor Supply','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(2,3,6);
hold on;
plot(vTime,wage,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Wage','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;
