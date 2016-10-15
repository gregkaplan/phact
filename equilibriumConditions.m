function vResidual = equilibriumConditions(vars)

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T ...
    nVars nEErrors rSS wSS KSS uSS cSS dSS VSS gSS dab varsSS

%% UNPACK VARIABLES

V = vars(1:I*J*N) + varsSS(1:I*J*N); % value function comes first

g = vars(I*J*N+1:2*I*J*N-1) + varsSS(I*J*N+1:2*I*J*N-1); % distribution comes second, last point is redundant
dab_aux = reshape(dab,I*J*N,1); dab_aux = dab_aux(1:end-1);
g_end = (1 - sum(g .* dab_aux))/dab(I,J,N);
gg = [g; g_end];

K = vars(2*I*J*N) + varsSS(2*I*J*N);
r_a = vars(2*I*J*N+1) + varsSS(2*I*J*N+1);
w = vars(2*I*J*N+2) + varsSS(2*I*J*N+2);
tfp = vars(2*I*J*N+3) + varsSS(2*I*J*N+3);

VDot = vars(nVars + 1: nVars + I*J*N);
gDot = vars(nVars + I*J*N + 1: nVars + 2*I*J*N-1);
tfpdot = vars(nVars + 2*I*J*N+3);

VEErrors = vars(2*nVars + 1: 2*nVars + I*J*N);
tfpshock = vars(2*nVars + I*J*N + 1);

KL = K;
V  = reshape(V,I,J,N);

%% ONE ITERATION OF HJB

% grids

set_grids
r_a_grid = repmat(r_a,I,J,N);
[l_grid,~,~] = derive_aggregates(KL);

% derivatives w.r.t. a
% preparations
VaF = 0*V;
VaB = 0*V;
%Vamin = myAD(0);
Vamin = 0;

% forward difference
VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))./daf_grid(:,1:J-1,:);
VaF(:,1:J-1,:) = max(VaF(:,1:J-1,:),Vamin);

% backward difference
VaB(:,2:J,:) = (V(:,2:J,:)-V(:,1:J-1,:))./dab_grid(:,2:J,:);
VaB(:,2:J,:) = max(VaB(:,2:J,:),Vamin);

% derivatives w.r.t. b
% preparations
VbF = 0*V;
VbB = 0*V;
Vbmin = 1e-8;

% forward difference
VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./dbf_grid(1:I-1,:,:);
VbF(1:I-1,:,:) = max(VbF(1:I-1,:,:),Vbmin);

% backward difference
VbB(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./dbb_grid(2:I,:,:);
VbB(2:I,:,:) = max(VbB(2:I,:,:),Vbmin);

% consumption decision  
% step 0: preparations
cF = 0 * V;
sF = 0 * V;
HcF = 0 * V;

cB = 0 * V;
sB = 0 * V;
HcB = 0 * V;

c0 = 0 * V;
s0 = 0 * V;
Hc0 = 0 * V;

% step 1: forward, backward and steady consumption and savings
cF(1:I-1,:,:) = VbF(1:I-1,:,:).^(-1/gamma);
cF(I,:,:) = zeros(1,J,N);
sF(1:I-1,:,:) = ((1-xi)-tau_I) * w * l_grid(1:I-1,:,:) .* y_grid(1:I-1,:,:) + b_grid(1:I-1,:,:) .* (r_b_grid(1:I-1,:,:) + deathrate*pam) + trans_grid(1:I-1,:,:) - cF(1:I-1,:,:);
sF(I,:,:) = zeros(1,J,N);
HcF(1:I-1,:,:) = u_fn(cF(1:I-1,:,:),gamma) + VbF(1:I-1,:,:) .* sF(1:I-1,:,:);
HcF(I,:,:) = -1e12*ones(1,J,N);
validF = (sF > 0);

cB(2:I,:,:) = VbB(2:I,:,:).^(-1/gamma);
cB(1,:,:) = ((1-xi)-tau_I) * w * l_grid(1,:,:) .* y_grid(1,:,:) + b_grid(1,:,:) .* (r_b_grid(1,:,:) + deathrate*pam) + trans_grid(1,:,:);
sB(2:I,:,:) = ((1-xi)-tau_I) * w * l_grid(2:I,:,:) .* y_grid(2:I,:,:) + b_grid(2:I,:,:) .* (r_b_grid(2:I,:,:) + deathrate*pam) + trans_grid(2:I,:,:) - cB(2:I,:,:);
sB(1,:,:) = zeros(1,J,N);
HcB(:,:,:) = u_fn(cB,gamma) + VbB .* sB;
validB = (sB < 0);

c0(:,:,:) = ((1-xi)-tau_I) * w * l_grid(:,:,:) .* y_grid(:,:,:) + b_grid(:,:,:) .* (r_b_grid(:,:,:) + deathrate*pam) + trans_grid(:,:,:);
s0(:,:,:) = zeros(I,J,N);
Hc0(:,:,:) = u_fn(c0,gamma);

% step 2: indicators of which version to use
IcF = validF .* max(~validB,(HcF>=HcB)) .* (HcF>=Hc0);
IcB = validB .* max(~validF,(HcB>=HcF)) .* (HcB>=Hc0);
Ic0 = 1 - IcF - IcB;

% step3: optimal consumption and savings
c = IcF .* cF + IcB .* cB + Ic0 .* c0;
s = IcF .* sF + IcB .* sB + Ic0 .* s0;

% deposit decision
% step 0: preparations
dFB = 0*V;
HdFB = 0*V;

dBF = 0*V;
HdBF = 0*V;

dBB = 0*V;
HdBB = 0*V;

% step 1: all possible mixture deposits
dFB(2:I,1:J-1,:) = opt_deposits(VaF(2:I,1:J-1,:),VbB(2:I,1:J-1,:),a_grid(2:I,1:J-1,:));
dFB(:,J,:) = zeros(I,1,N);
dFB(1,1:J-1,:) = zeros(1,J-1,N);
HdFB(2:I,1:J-1,:) = VaF(2:I,1:J-1,:) .* dFB(2:I,1:J-1,:) - VbB(2:I,1:J-1,:) .* (dFB(2:I,1:J-1,:) + adj_cost_fn(dFB(2:I,1:J-1,:),a_grid(2:I,1:J-1,:)));
HdFB(:,J,:) = -1.0e12 * ones(I,1,N);
HdFB(1,1:J-1,:) = -1.0e12 * ones(1,J-1,N);
validFB = (dFB > 0) .* (HdFB > 0);

dBF(1:I-1,2:J,:) = opt_deposits(VaB(1:I-1,2:J,:),VbF(1:I-1,2:J,:),a_grid(1:I-1,2:J,:));
dBF(:,1,:) = zeros(I,1,N);
dBF(I,2:J,:) = zeros(1,J-1,N);
HdBF(1:I-1,2:J,:) = VaB(1:I-1,2:J,:) .* dBF(1:I-1,2:J,:) - VbF(1:I-1,2:J,:) .* (dBF(1:I-1,2:J,:) + adj_cost_fn(dBF(1:I-1,2:J,:),a_grid(1:I-1,2:J,:)));
HdBF(:,1,:) = -1.0e12 * ones(I,1,N);
HdBF(I,2:J,:) = -1.0e12 * ones(1,J-1,N);
validBF = (dBF <= -adj_cost_fn(dBF,a_grid)) .* (HdBF > 0);

VbB(1,2:J,:) = u_fn(cB(1,2:J,:),gamma);
dBB(:,2:J,:) = opt_deposits(VaB(:,2:J,:),VbB(:,2:J,:),a_grid(:,2:J,:));
dBB(:,1,:) = zeros(I,1,N);
HdBB(:,2:J,:) = VaB(:,2:J,:) .* dBB(:,2:J,:) - VbB(:,2:J,:) .* (dBB(:,2:J,:) + adj_cost_fn(dBB(:,2:J,:),a_grid(:,2:J,:)));
HdBB(:,1,:) = -1.0e12 * ones(I,1,N);
validBB = (dBB > -adj_cost_fn(dBB,a_grid)) .* (dBB <= 0) .* (HdBB > 0);

% step 2: indicators of which version to use
IcFB = validFB .* max(~validBF,(HdFB >= HdBF)) .* max(~validBB,(HdFB >= HdBB));
IcBF = max(~validFB,(HdBF >= HdFB)) .* validBF .* max(~validBB,(HdBF >= HdBB));
IcBB = max(~validFB,(HdBB >= HdFB)) .* max(~validBF,(HdBB >= HdBF)) .* validBB;
Ic00 = (~validFB) .* (~validBF) .* (~validBB);

% step 3: optimal deposits
d = IcFB .* dFB + IcBF .* dBF + IcBB .* dBB + Ic00 .* zeros(I,J,N);

% compute various objects of interest
u = u_fn(c,gamma);

% construct transition matrix
% preparations   
% placeholders

audriftB = 0*V;
audriftF = 0*V;
budriftB = 0*V;
budriftF = 0*V;

chi = 0*V;
yy = 0*V;
zeta = 0*V;

X = 0*V;
Y = 0*V;
Z = 0*V;

chiu = 0*V;
yyu = 0*V;
zetau = 0*V;

Xu = 0*V;
Yu = 0*V;
Zu = 0*V;

% compute all drifts
adriftB = min(d,0) + min(a_grid .* (r_a_grid + deathrate*pam) + xi * w * l_grid .* y_grid,0);
adriftF = max(d,0) + max(a_grid .* (r_a_grid + deathrate*pam) + xi * w * l_grid .* y_grid,0);

bdriftB = min(-d - adj_cost_fn(d,a_grid),0) + min(s,0);
bdriftF = max(-d - adj_cost_fn(d,a_grid),0) + max(s,0);

audriftB(1:I-1,:,:) = min(d(1:I-1,:,:) + a_grid(1:I-1,:,:) .* (r_a_grid(1:I-1,:,:) + deathrate*pam) + xi * w * l_grid(1:I-1,:,:) .* y_grid(1:I-1,:,:),0);
audriftB(I,:,:) = min(d(I,:,:) + a_grid(I,:,:) .* (r_a_grid(I,:,:) + deathrate*pam) + xi * w * l_grid(I,:,:) .* y_grid(I,:,:),0);
audriftF(1:I-1,:,:) = max(d(1:I-1,:,:) + a_grid(1:I-1,:,:) .* (r_a_grid(1:I-1,:,:) + deathrate*pam) + xi * w * l_grid(1:I-1,:,:) .* y_grid(1:I-1,:,:),0);
audriftF(I,:,:) = max(d(I,:,:) + a_grid(I,:,:) .* (r_a_grid(I,:,:) + deathrate*pam) + xi * w * l_grid(I,:,:) .* y_grid(I,:,:),0);

budriftB(1:I-1,:,:) = min(s(1:I-1,:,:) - d(1:I-1,:,:) - adj_cost_fn(d(1:I-1,:,:),a_grid(1:I-1,:,:)),0);
budriftB(I,:,:) = min(s(I,:,:) - d(I,:,:) - adj_cost_fn(d(I,:,:),a_grid(I,:,:)),0);
budriftF(1:I-1,:,:) = max(s(1:I-1,:,:) - d(1:I-1,:,:) - adj_cost_fn(d(1:I-1,:,:),a_grid(1:I-1,:,:)),0);
budriftF(I,:,:) = max(s(I,:,:) - d(I,:,:) - adj_cost_fn(d(I,:,:),a_grid(I,:,:)),0);   

% transition matrix for a
chi(:,2:J,:) = -adriftB(:,2:J,:)./dab_grid(:,2:J,:); chi(:,1,:) = zeros(I,1,N);
yy(:,2:J-1,:) = adriftB(:,2:J-1,:)./dab_grid(:,2:J-1,:) - adriftF(:,2:J-1,:)./daf_grid(:,2:J-1,:); yy(:,1,:) = - adriftF(:,1,:)./daf_grid(:,1,:); yy(:,J,:) = adriftB(:,J,:)./dab_grid(:,J,:);
zeta(:,1:J-1,:) = adriftF(:,1:J-1,:)./daf_grid(:,1:J-1,:); zeta(:,J,:) = zeros(I,1,N);

centdiag = reshape(yy,I*J,N);
lowdiag = reshape(chi,I*J,N);
lowdiag = circshift(lowdiag,-I);
updiag = reshape(zeta,I*J,N);
updiag = circshift(updiag,I);

centdiag = reshape(centdiag,I*J*N,1);
updiag   = reshape(updiag,I*J*N,1);
lowdiag  = reshape(lowdiag,I*J*N,1);

aa = spdiags(centdiag,0,I*J*N,I*J*N) + spdiags(updiag,I,I*J*N,I*J*N) + spdiags(lowdiag,-I,I*J*N,I*J*N);

chiu(:,2:J,:) = -audriftB(:,2:J,:)./dab_grid(:,2:J,:); chiu(:,1,:) = zeros(I,1,N);
yyu(:,2:J-1,:) = audriftB(:,2:J-1,:)./dab_grid(:,2:J-1,:) - audriftF(:,2:J-1,:)./daf_grid(:,2:J-1,:); yyu(:,1,:) = - audriftF(:,1,:)./daf_grid(:,1,:); yyu(:,J,:) = audriftB(:,J,:)./dab_grid(:,J,:);
zetau(:,1:J-1,:) = audriftF(:,1:J-1,:)./daf_grid(:,1:J-1,:); zetau(:,J,:) = zeros(I,1,N);

centdiagu = reshape(yyu,I*J,N);
lowdiagu = reshape(chiu,I*J,N);
lowdiagu = circshift(lowdiagu,-I);
updiagu = reshape(zetau,I*J,N);
updiagu = circshift(updiagu,I);

centdiagu = reshape(centdiagu,I*J*N,1);
updiagu   = reshape(updiagu,I*J*N,1);
lowdiagu  = reshape(lowdiagu,I*J*N,1);

aau = spdiags(centdiagu,0,I*J*N,I*J*N) + spdiags(updiagu,I,I*J*N,I*J*N) + spdiags(lowdiagu,-I,I*J*N,I*J*N);
centdiag_aa = centdiagu;
updiag_aa = updiagu;
lowdiag_aa = lowdiagu;

% transition matrix for b
X(2:I,:,:) = - bdriftB(2:I,:,:)./dbb_grid(2:I,:,:); X(1,:,:) = zeros(1,J,N);
Y(2:I-1,:,:) = bdriftB(2:I-1,:,:)./dbb_grid(2:I-1,:,:) - bdriftF(2:I-1,:,:)./dbf_grid(2:I-1,:,:); Y(1,:,:) = - bdriftF(1,:,:)./dbf_grid(1,:,:); Y(I,:,:) = bdriftB(I,:,:)./dbb_grid(I,:,:);
Z(1:I-1,:,:) = bdriftF(1:I-1,:,:)./dbf_grid(1:I-1,:,:); Z(I,:,:) = zeros(1,J,N);

centdiag = reshape(Y,I*J,N);
lowdiag = reshape(X,I*J,N);
lowdiag = circshift(lowdiag,-1);
updiag = reshape(Z,I*J,N);
updiag = circshift(updiag,1);

centdiag = reshape(centdiag,I*J*N,1);
updiag   = reshape(updiag,I*J*N,1);
lowdiag  = reshape(lowdiag,I*J*N,1);

bb = spdiags(centdiag,0,I*J*N,I*J*N) + spdiags(updiag,1,I*J*N,I*J*N) + spdiags(lowdiag,-1,I*J*N,I*J*N);

Xu(2:I,:,:) = - budriftB(2:I,:,:)./dbb_grid(2:I,:,:); Xu(1,:,:) = zeros(1,J,N);
Yu(2:I-1,:,:) = budriftB(2:I-1,:,:)./dbb_grid(2:I-1,:,:) - budriftF(2:I-1,:,:)./dbf_grid(2:I-1,:,:); Yu(1,:,:) = - budriftF(1,:,:)./dbf_grid(1,:,:); Yu(I,:,:) = budriftB(I,:,:)./dbb_grid(I,:,:);
Zu(1:I-1,:,:) = budriftF(1:I-1,:,:)./dbf_grid(1:I-1,:,:); Zu(I,:,:) = zeros(1,J,N);

centdiagu = reshape(Yu,I*J,N);
lowdiagu = reshape(Xu,I*J,N);
lowdiagu = circshift(lowdiagu,-1);
updiagu = reshape(Zu,I*J,N);
updiagu = circshift(updiagu,1);

centdiagu = reshape(centdiagu,I*J*N,1);
updiagu   = reshape(updiagu,I*J*N,1);
lowdiagu  = reshape(lowdiagu,I*J*N,1);

bbu = spdiags(centdiagu,0,I*J*N,I*J*N) + spdiags(updiagu,1,I*J*N,I*J*N) + spdiags(lowdiagu,-1,I*J*N,I*J*N);

centdiag_bb = centdiagu;
updiag_bb = updiagu;
lowdiag_bb = lowdiagu;

% transition matrix for productivity state
cc = kron(lambda,speye(I*J));
ccu = kron(lambda,speye(I*J));

% full transition matrix

A = aa + bb + cc;

%% COMPUTE EQUILIBRIUM CONDITIONS

% HJB Equation
hjbResidual = reshape(u,I*J*N,1) + A * reshape(V,I*J*N,1) + VDot + VEErrors - (rho + deathrate) * reshape(V,I*J*N,1);

% KFE
A = aau + bbu + ccu;
A = A';
gg_tilde = dab_tilde_mat * gg;
gIntermediate = A * gg_tilde;
dab_tilde_mat_inv = spdiags(repmat(1./dab_tilde,N,1),0,I*J*N,I*J*N);
% gIntermediate = dab_tilde_mat \ gIntermediate;
gIntermediate = dab_tilde_mat_inv * gIntermediate;

dab_small = reshape(dab(:,:,1),I*J,1);

loc = find(b==0);
dab_small = dab_small./dab_small(loc)*deathrate;
dab_small(loc) = 0;
death_process = -deathrate*speye(I*J);
death_process(loc,:)=dab_small;
death_process = kron(speye(N),death_process);

gIntermediate = gIntermediate + death_process * gg;

gResidual = gDot - gIntermediate(1:I*J*N-1,1);

% Aggregate conditions
K_upd = sum(reshape(a_grid,I*J*N,1) .* gg .* reshape(dab,I*J*N,1));
kResidual = K_upd - K;
raResidual = exp(tfp) * alpha * (K_upd^(alpha - 1)) - delta - r_a;
wResidual = exp(tfp) * (1 - alpha) * (K_upd^alpha) - w;

% Law of motion for aggregate shocks
tfpResidual = tfpdot + (1 - rho_tfp) * tfp - sigma_tfp * tfpshock;

vResidual = [hjbResidual;gResidual;kResidual;raResidual;wResidual;tfpResidual];