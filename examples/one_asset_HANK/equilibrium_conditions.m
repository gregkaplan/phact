% Unpack vars
V  = reshape(vars(1:I*J,1) + vars_SS(1:I*J,1),I,J);
inflation = vars(n_v) + vars_SS(n_v);

gg = vars(n_v+1:n_v+n_g-1) + vars_SS(n_v+1:n_v+n_g-1);
azdeltavec = reshape(azdelta,I*J,1);
g_end = (1 - sum(gg .* azdeltavec(1:end-1))) ./ azdeltavec(end) ;
g = [gg;g_end];
MP = vars(n_v+n_g,1) + vars_SS(n_v+n_g);

w = vars(n_v+n_g+1) + vars_SS(n_v+n_g+1);
hours = vars(n_v+n_g+2) + vars_SS(n_v+n_g+2);
consumption = vars(n_v+n_g+3) + vars_SS(n_v+n_g+3);
output = vars(n_v+n_g+4) + vars_SS(n_v+n_g+4);
assets = vars(n_v+n_g+5) + vars_SS(n_v+n_g+5);

VDot = vars(nVars+1:nVars+n_v-1);
inflationDot = vars(nVars+n_v);
gDot = vars(nVars+n_v+1:nVars+n_v+n_g-1);
MPDot = vars(nVars+n_v+n_g);

VEErrors = vars((2*nVars+1):(2*nVars+n_v-1),1);
inflationEError = vars(2*nVars+n_v,1);

MPShock = vars(2*nVars+n_v+1,1);

TFP = 1;

% Initialize other variables, using V to ensure everything is a dual number
Vaf = V;
Vab = V;

m = w/TFP;
profit = (1 - m) .* output;
profshare = zz./meanlabeff .* profit;
r_Nominal = r_SS + taylor_inflation * inflation + taylor_outputgap * (log(output) - log(Y_SS)) + MP;
r = r_Nominal - inflation;
lumptransfer = labtax * w * hours - G_SS - (r_Nominal - (1-govbcrule_fixnomB) * inflation) * assets;

%% HJB
h0 = h_SS;
h = h_SS;

% forward difference
Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:)) ./ dazf(1:I-1,:);
Vaf(I,:) = (h(I,:) .* zz(I,:) .* w .* (1-labtax) + lumptransfer + profshare(I,:) + r.*amax ).^(-coefrra);

% backward difference
Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:)) ./ dazb(2:I,:);
Vab(1,:) = (h0(1,:) .* zz(1,:) .* w .* (1-labtax) + lumptransfer + profshare(1,:) + r.*amin ).^(-coefrra);

% consumption and savings with forward difference
cf = Vaf.^(-1/coefrra);
hf = (zz .* w .* (1-labtax) .* Vaf ./ labdisutil ) .^ frisch;
hf = min(hf,maxhours);

for ih = 1:niter_hours
	cf(I,:)  = hf(I,:) .* zz(I,:) .* w .* (1-labtax) + lumptransfer + profshare(I,:) + r.*aa(I,:);
	hf(I,:) = (zz(I,:) .* w .* (1-labtax) .* (cf(I,:).^(-coefrra))./ labdisutil ) .^ frisch;
	hf(I,:) = min(hf(I,:),maxhours);
end
Vaf(I,:) = cf(I,:).^(-coefrra);

sf = hf .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa - cf;
if coefrra==1
	Vf = log(cf) - labdisutil .* (hf.^(1+1/frisch)/(1+1/frisch));
else
	Vf = (cf.^(1-coefrra))./(1-coefrra) - labdisutil .* (hf.^(1+1/frisch)/(1+1/frisch));
end
Vf = (cf>0).*(Vf + sf.*Vaf) + (cf<=0).*(-1e12);

% consumption and savings with backward difference
cb = Vab.^(-1/coefrra);
hb = (zz .* w .* (1-labtax) .* Vab ./ labdisutil ) .^ frisch;
hb = min(hb,maxhours);

for ih = 1:niter_hours
	cb(1,:)  = hb(1,:) .* zz(1,:) .* w .* (1-labtax) + lumptransfer + profshare(1,:) + r.*aa(1,:);
	hb(1,:) = (zz(1,:) .* w .* (1-labtax) .* (cb(1,:).^(-coefrra))./ labdisutil ) .^ frisch;
	hb(1,:) = min(hb(1,:),maxhours);
end
Vab(1,:) = cb(1,:).^(-coefrra);

sb = hb .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa - cb;
if coefrra==1
	Vb = log(cb) - labdisutil .* (hb.^(1+1/frisch)/(1+1/frisch));
else
	Vb = (cb.^(1-coefrra))./(1-coefrra) - labdisutil .* (hb.^(1+1/frisch)/(1+1/frisch));
end
Vb = (cb>0).*(Vb + sb.*Vab) + (cb<=0).*(-1e12);

% consumption and derivative of value function at steady state
for ih = 1:niter_hours
	c0 = h0 .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa;
	h0 = (zz .* w .* (1-labtax) .* (c0.^(-coefrra)) ./ labdisutil ) .^ frisch;
	h0 = min(h0,maxhours);
end
c0 = h0 .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa;

Va0 = (c0).^(-coefrra);
if coefrra==1
	V0 = log(c0) - labdisutil .* (h0.^(1+1/frisch)/(1+1/frisch));
else
	V0 = (c0.^(1-coefrra))./(1-coefrra) - labdisutil .* (h0.^(1+1/frisch)/(1+1/frisch));
end
V0 = (c0>0).*V0 + (c0<=0).*(-1e12);

% choose upwinding direction
Ineither = (1-(sf>0)) .* (1-(sb<0));
Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
Iboth = (sb<0).*(sf>0);
Ib = Iunique.*(sb<0).*(Vb>V0) + Iboth.*(Vb==max(max(Vb,Vf),V0));
If = Iunique.*(sf>0).*(Vf>V0) + Iboth.*(Vf==max(max(Vb,Vf),V0));
I0 = Ineither + (1-Ineither).*(V0==max(max(Vb,Vf),V0));
I0 = 1-Ib-If;

Va = Vaf.*If + Vab.*Ib + Va0.*I0;
h = hf.*If + hb.*Ib + h0.*I0;
c = cf.*If + cb.*Ib + c0.*I0;
s = sf.*If + sb.*Ib;

if coefrra==1
	u = log(c) - labdisutil * (h.^(1 + 1/frisch))/(1+1/frisch);
else
	u = (c.^(1-coefrra))./(1-coefrra) - labdisutil .* (h.^(1 + 1/frisch))/(1+1/frisch);
end

% CONSTRUCT MATRIX A
X = -Ib.*sb./dazb;
Y = -If.*sf./dazf + Ib.*sb./dazb;
Z = If.*sf./dazf;

updiag=0; %This is needed because of the peculiarity of spdiags.
for j=1:J
    updiag=[updiag;Z(1:I-1,j);0];
end

centdiag=reshape(Y,I*J,1);

lowdiag=X(2:I,1);
for j=2:J
    lowdiag=[lowdiag;0;X(2:I,j)];
end

AA = spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);

A = AA + Aswitch;

%% Collect dynamics
% HJB equation
hjbResidual = reshape(u,I*J,1) + A * reshape(V,I*J,1) + VDot + VEErrors - rrho * reshape(V,I*J,1);

% Inflation
pcResidual = -((r - 0) * inflation - (ceselast/priceadjust * (w/TFP - (ceselast-1)/ceselast) + inflationDot - inflationEError));

% KFE
gIntermediate = spdiags(1./diag(azdelta_mat),0,I*J,I*J) * A' * (azdelta_mat * g);
gResidual = gDot - gIntermediate(1:end-1,1);

% Monetary Policy
MPResidual = MPDot - (-ttheta_MP * MP + ssigma_MP * MPShock);

% Market Clearing
realsav = sum(aa(:).*g(:).*azdelta(:));
realsavDot = sum(s(:).*g(:).*azdelta(:));
bondmarketResidual = realsavDot/realsav + govbcrule_fixnomB * inflation;
labmarketResidual = sum(zz(:).* h(:) .* g(:) .* azdelta(:)) - hours;
consumptionResidual = sum(c(:).*g(:).*azdelta(:)) - consumption;
outputResidual = TFP * hours - output;
assetsResidual = assets - realsav;

v_residual = [hjbResidual;
              pcResidual;
              gResidual;
              MPResidual;
              bondmarketResidual;
              labmarketResidual;
              consumptionResidual;
              outputResidual;
              assetsResidual];
