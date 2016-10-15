% Inputs:  (1) vars: vector which contains the variables in the system, the time derivatives of 
%			         those variables, the expectational errors, and the shocks
%
% Outputs: (1) vResduals: residuals of equilibrium conditions, evaluated at vars

function vResidual = equilibriumConditions_flex(vars)

%----------------------------------------------------------------
% Housekeeping
%----------------------------------------------------------------

% Declare global variables
global ggamma rrho pphi ddeath cchi nSS rrho0 rrhoMin rrhoMax wedge eepsilon ttheta pphiInflation pphiOutput nnuTFP ...
	ssigmaTFP nnuMon ssigmaMon z J z_avg g_z I I_neg zz zzz Aswitch rmin rmax r0 maxit crit Delta Ir crit_S ...
	Irrho crit_objective T dt N nVars nEErrors nShocks aalpha ddelta amin amax a aa daa_f daa_b nVars YSS ...
	IfSS IbSS I0SS rSS da_stacked da_tilde varsSS ttau ttransfer
	
% Unpack vars
V = reshape(vars(1:I*J,1) + varsSS(1:I*J,1),I,J);
gg = vars(I*J+1:2*I*J-1) + varsSS(I*J+1:2*I*J-1);
g_end = (1 / da_stacked(I*J,1)) * (1 - sum(gg .* da_stacked(1:I*J-1)));
g = [gg;g_end];
logAggregateTFP = vars(2*I*J,1);

VDot = vars(nVars+1:nVars+I*J,1);
gDot = vars(nVars+I*J+1:nVars+2*I*J-1,1);
logAggregateTFPDot = vars(nVars+2*I*J,1);

VEErrors = vars(2*nVars+1:2*nVars+I*J,1);

logAggregateTFPShock = vars(2*nVars+nEErrors+1,1);

% Initialize other variables, using V to ensure everything is a dual number
dVf = V;
dVb = V;

%----------------------------------------------------------------
% Compute auxiliary variables from known relationships
%----------------------------------------------------------------

K = (reshape(aa,I*J,1) .* g)'*da_stacked;
N = (((1 - ttau) / cchi) * (z_avg ^ pphi) * (1 - aalpha) * exp(logAggregateTFP) * (K ^ aalpha)) ^ (1 / (pphi + aalpha));
output = exp(logAggregateTFP) * (K ^ aalpha) * (N ^ (1 - aalpha));
r_real = aalpha * output / K - ddelta;
w = (1 - aalpha) * output / N;
nHours = N / z_avg;

% Adjust rates of return by death
rrho_effective = rrho + ddeath;
r = ddeath + r_real .* (aa >= 0) + (r_real + wedge) .* (aa < 0);

%----------------------------------------------------------------
% HJB Equation
%----------------------------------------------------------------

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

% Construct matrix for implicit updating scheme
X = -sb .* IbSS ./ daa_b;
Y = -sf .* IfSS ./ daa_f + sb .* IbSS ./ daa_b;
Z = sf .* IfSS ./ daa_f;

% The following is needed because of a peculiarity of spdiags
updiag = [0];
for j=1:J
	updiag=[updiag;Z(1:I-1,j);0];
end
centdiag = reshape(Y,I*J,1);
lowdiag=X(2:I,1);
for j=2:J
	lowdiag=[lowdiag;0;X(2:I,j)];
end

% Create A matrix
AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
A = AA + Aswitch;

%----------------------------------------------------------------
% KFE
%----------------------------------------------------------------

% Mass preservation
tmp = repmat(da_stacked,1,I*J);
A_adj = (tmp .* A) ./ (tmp');

% Adjust for death
aux = sparse(I,I);
aux(:,1)=ddeath*da_tilde/da_tilde(1);
aux2 = kron(speye(J),aux);
A_adj = A_adj + aux2 - ddeath*speye(I*J);

% Updating matrix
AT = A_adj';

%----------------------------------------------------------------
% Compute equilibrium conditions
%----------------------------------------------------------------

% HJB equation
hjbResidual = reshape(u,I*J,1) + A * reshape(V,I*J,1) + VDot - VEErrors - rrho_effective * reshape(V,I*J,1);

% KFE
gDotIntermediate = AT * g;
gResidual = gDot - gDotIntermediate(1:I*J-1,1);

% Laws of motion for aggregate shocks
tfpResidual = logAggregateTFPDot - (-nnuTFP * logAggregateTFP + ssigmaTFP * logAggregateTFPShock);

% Package
vResidual = [hjbResidual;gResidual;tfpResidual];
