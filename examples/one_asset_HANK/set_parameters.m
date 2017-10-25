%% Economic Parameters

% Preferences
coefrra = 1;
frisch = 0.5;
adjfricshgridfrac = 0.85;
meanlabeff = 3; 	%so that at h=1/3 outout will be approximately = 1;
labdisutil = meanlabeff ./ ( (0.75 .^(-coefrra)) * ((1.0/3.0).^(1.0/frisch))); %guess labdisutil so that at average wages and average consumption hours =1/3 (sets C/Y = 0.75);
maxhours = 1;

% Production
ceselast = 10;							 % elasticity of substitution / demand
priceadjust = 100;

% Policy parameters
taylor_inflation = 1.25;	%taylor rule coefficient on inflation
taylor_outputgap = 0;		%taylor rule coefficient on output
labtax = 0.2;			% marginal tax rate
govbondtarget = 6;		%multiple of quarterly GDP
lumptransferpc = 0.06;	 %6% of quarterly GDP in steady state
govbcrule_fixnomB = 0.0;

% Some aggregates
Y_SS = 1;
B_SS = govbondtarget .* Y_SS;
m_SS = (ceselast-1)/ceselast;
w_SS = m_SS;
N_SS = 1/3;						    % steady state hours: so that quartelry GDP = 1 in s.s
profit_SS = (1 - m_SS) * Y_SS;
lumptransfer = lumptransferpc .* Y_SS;

% Aggregate shocks
ssigma_MP = sqrt(0.05);
ttheta_MP = 0.25;


%% Approximation Parameters
% HJB
maxit_HJB = 500;
tol_HJB  = 1.0e-8;
Delta_HJB = 1.0e6;

maxit_KFE = 1000;
tol_KFE  = 1.0e-12;
Delta_KFE = 1.0e6;

crit  = 10^(-6);

% steady state r
r0 = 0.005;
rmax = 0.08;
rmin = 0.001;
Ir = 100;

% steady state rho iteration
rrho0 = 0.02;
rhomax = 0.05;
rhomin = 0.005;
crit_S = 10^(-5);

% hours loop
niter_hours = 10;


%% Set Grids

%% Asset
I          = 100;
agridparam = 1; % 1 for linear
amin       = 0;
amax       = 40;

a          = linspace(0,1,I)';
a          = a.^(1/agridparam);
a          = amin + (amax-amin).*a;

%% Income
J = 2;

ygrid_combined = [0.2; 1];
ymarkov_combined = [-0.5 0.5; 0.0376 -0.0376];

% Compute stationary income distribution
AT = ymarkov_combined';
g_z = ones(J,1)/J;
for n = 1 : 50
	g_z_new = (speye(J) - AT * 1000) \ g_z;
	diff = max(abs(g_z_new - g_z));
	if diff < 10e-6
		break
	end
	g_z = g_z_new;
end

z = exp(ygrid_combined);

% scale so that mean = meanlabeff
z  	  = meanlabeff.*z./sum(z .* g_z);
z_bar = sum(z .* g_z);
z     = z';
zz = ones(I,1)*z;


%% Kronecker to matrix
aa = a*ones(1,J);

daf = [a(2:I)-a(1:I-1); a(I)-a(I-1)];
dazf = repmat(daf,1,J);
dab = [a(2)-a(1); a(2:I)-a(1:I-1)];
dazb = repmat(dab,1,J);

adelta = zeros(size(a));
adelta(1) = 0.5.*daf(1);
adelta(2:I-1) = 0.5.*(daf(1:I-2)+daf(2:I-1));
adelta(I) = 0.5*daf(I-1);

azdelta = repmat(adelta,1,J);
azdelta_mat = spdiags(azdelta(:),0,I*J,I*J);

Aswitch = kron(ymarkov_combined,speye(I));


%% Others
IterateR = 0;
IterateRho = 1;
