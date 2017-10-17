%% Solves the (non-stochastic) steady state of one-asset HANK model

% Initialize matrices for finite differences
Vaf = zeros(I,J);
Vab = zeros(I,J);

% Initialize interest rate
r = r0;
rrho = rrho0;

profshare = zz./meanlabeff * profit_SS;
lumptransfer = lumptransferpc .* Y_SS;

for ir = 1:Ir
	c = zeros(I,J);
	h = ones(I,J).*(1/3);
	h0 = ones(I,J);

	w = w_SS;

	% initial guess
	if coefrra==1
		v0 = (log(h .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa ) - labdisutil .* h.^(1+1/frisch)/(1+1/frisch) )./rrho;
	else
		v0 = ((h .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa ).^(1-coefrra)./(1-coefrra) - labdisutil .* (h.^(1+1/frisch))./(1+1/frisch) ) ./rrho;
	end
	v = v0;

	% iterate HJB
	for ihjb = 1:maxit_HJB
	    V = v;

	    % forward difference
	    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:)) ./ dazf(1:I-1,:);
	    Vaf(I,:) = (h(I,:) .* zz(I,:) .* w .* (1-labtax) + lumptransfer + profshare(I,:) + r.*amax ).^(-coefrra);

	    % backward difference
	    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))  ./ dazb(2:I,:);
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

	    AA =spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);

	    A = AA + Aswitch;
	    B = (1/Delta_HJB + rrho)*speye(I*J) - A;

	    u_stacked = reshape(u,I*J,1);
	    V_stacked = reshape(V,I*J,1);

	    b = u_stacked + V_stacked/Delta_HJB;

	    V_stacked = B\b;

	    V = reshape(V_stacked,I,J);

	    Vchange = V - v;
	    v = V;

	    err_HJB = max(max(abs(Vchange)));
		if DisplayLev >= 2
			disp(['HJB iteration ' int2str(ihjb) ', error = ' num2str(err_HJB)]);
		end
	    if err_HJB<tol_HJB
	        break
	    end
	end

	% solve KFE

	% initialize KFE iterations at a=0;
	g0 = zeros(I,J);
	g0(a==0,:) = g_z;
	g0 = g0./azdelta;

	gg = reshape(g0,I*J,1);

	% Solve linear system
	for ikfe = 1:maxit_KFE
		gg_tilde = azdelta_mat * gg;
		gg1_tilde = (speye(I*J) - Delta_KFE * A') \ gg_tilde;
		gg1_tilde = gg1_tilde ./ sum(gg1_tilde);
		gg1 = azdelta_mat\gg1_tilde;

		err_KFE = max(abs(gg1-gg));
		if DisplayLev >= 2
 			disp(['KFE iteration ' int2str(ikfe) ', distance is ' num2str(err_KFE)]);
		end
		if err_KFE < tol_KFE
 			if DisplayLev >= 2
				disp(['Distribution Converged, Iteration = ' int2str(ikfe)]);
			end
			break
		end

		gg = gg1;
	end

	g = reshape(gg,I,J);

	% equilibrium objects
	N_SS = sum(zz(:).* h(:) .* g(:) .* azdelta(:));
	Y_SS = exp(0) * N_SS ;
	B_SS = sum(g(:) .* aa(:) .* azdelta(:));
	m_SS = (ceselast-1)/ceselast;
	w_SS = m_SS;
	profit_SS = (1 - m_SS) .* Y_SS;
	profshare = zz./meanlabeff .* profit_SS;
	lumptransfer = lumptransferpc .* Y_SS;
	bond_err = B_SS./Y_SS - govbondtarget;

    if bond_err>crit_S
		if IterateR==1
			if DisplayLev >= 1
				disp(['Interest rate: ' num2str(r) ', Excess Supply: ' num2str(bond_err) ]);
			end
			rmax = r;
			r = 0.5*(r+rmin);
		elseif IterateRho==1
			if DisplayLev >= 1
				disp(['Discount rate: ' num2str(rrho) ', Excess Supply: ' num2str(bond_err) ]);
			end
			rhomin = rrho;
			rrho = 0.5*(rrho+rhomax);
		end

	elseif bond_err<-crit_S
		if IterateR==1
			if DisplayLev >= 1
				disp(['Interest rate: ' num2str(r) ', Excess Demand: ' num2str(bond_err) ]);
			end
			rmin = r;
			r = 0.5*(r+rmax);
		elseif IterateRho==1
			if DisplayLev >= 1
				disp(['Discount rate: ' num2str(rrho) ', Excess Demand: ' num2str(bond_err) ]);
			end
			rhomax = rrho;
			rrho = 0.5*(rrho+rhomin);
		end
    elseif abs(bond_err)<crit_S
			disp(['Steady State Found, Interest rate = ' num2str(r) ', Discount rate = ' num2str(rrho)]);
        break

    end

end

V_SS = V;
inflation_SS = 0;
g_SS = g;
r_SS = r;
u_SS = u;
c_SS = c;
h_SS = h;
s_SS = s;
rnom_SS = r_SS + inflation_SS;
B_SS = sum(g_SS(:) .* aa(:) .* azdelta(:));
N_SS = sum(zz(:).* h_SS(:) .* g_SS(:) .* azdelta(:));
Y_SS = exp(0) * N_SS ;
m_SS = (ceselast-1)/ceselast;
w_SS = m_SS;
profit_SS = (1 - m_SS) .* Y_SS;
C_SS = sum(c_SS(:) .* g_SS(:) .* azdelta(:));
T_SS = lumptransfer;
G_SS = labtax * w_SS * N_SS - T_SS - r_SS * B_SS;

% Collect savings and income (useful for checks)
% savings = h .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa - c;
% income = h .* zz .* w .* (1-labtax) + lumptransfer + profshare + r.*aa;

% Collect variables into vector
vars_SS = zeros(n_v + n_g + n_p, 1);

% Value Function
vars_SS(1:(n_v-1)) = reshape(V_SS,I*J, 1);

% Inflation
vars_SS(n_v,1) = inflation_SS;

% Distribution
gg_SS = reshape(g_SS,I*J,1);
vars_SS(n_v+1:n_v+n_g-1, 1) = gg_SS(1:I*J-1);
vars_SS(n_v+n_g, 1) = 0;
vars_SS(n_v+n_g+1, 1) = w_SS;
vars_SS(n_v+n_g+2, 1) = N_SS;
vars_SS(n_v+n_g+3, 1) = C_SS;
vars_SS(n_v+n_g+4, 1) = Y_SS;
vars_SS(n_v+n_g+5, 1) = B_SS;
