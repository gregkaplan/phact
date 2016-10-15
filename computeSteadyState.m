% Computes the steady state
%
% Outputs:  (1) r: steady state interest rate
%           (2) w: steady state wage
%			(3) rK: steady state capital rental rate
%			(4) Y: steady state output
%			(5) K: steady state capital stock
%			(6) m: steady state marginal cost
%           (7) V: value function over grid
%			(8) g: distribution
%           (9) dV_Upwind: derivative of value function by upwind scheme
%           (10) dVf: derivative of value function by forward difference
%           (11) dVb: derivative of value function by backward difference
%           (12) If: indicator for forward drift in savings
%           (13) Ib: indicator for backward drift in savings
%           (14) I0: indicator for no drift in savings

function [r_guess,w,rK,output,A_st,m,V,g,c,dV_Upwind,dVf,dVb,If,Ib,I0] = computeSteadyState()

%----------------------------------------------------------------
% Housekeeping
%----------------------------------------------------------------

% Declare global variables
global ggamma rrho pphi ddeath cchi nSS rrho0 rrhoMin rrhoMax wedge eepsilon ttheta pphiMon pphiOutput nnuTFP ...
	ssigmaTFP nnuMon ssigmaMon z J z_avg g_z I I_neg zz zzz Aswitch rmin rmax r0 maxit crit Delta Ir crit_S ...
	Irrho crit_objective T dt N nVars nEErrors nShocks aalpha ddelta amin amax a aa daa_f daa_b da_stacked ...
	da_tilde ttau z_cutoff ttransfer A
	
% Marginal cost in steaady state
m = (eepsilon - 1) / eepsilon;

% Initialize matrices for finite differences
dVbf = zeros(I,J);
dVbb = zeros(I,J);
dVzf = zeros(I,J);
dVzb = zeros(I,J);
dVzz = zeros(I,J);
c = zeros(I,J);

% Initialize discount rate for calibration
rrho = rrho0;
r_guess = r0;

%----------------------------------------------------------------
% Iterate to find discount rate which matches wealth to income ratio
%----------------------------------------------------------------

for iRrho = 1 : Irrho

	rrho_effective = rrho + ddeath;
	rmin = .0001;								% lower bound for steady state interest rate
	rmax = rrho_effective;						% upper bound for steady state interest rate

	for ir = 1 : Ir
	
		r_r(ir)=r_guess;
		rmin_r(ir)=rmin;
		rmax_r(ir)=rmax;
	
		%----------------------------------------------------------------
		% Preliminaries given real interest rate r
		%----------------------------------------------------------------
	
		% Given r, compute aggregates and other prices
		K = (((1 - (1 - aalpha) * m) / (r_guess + ddelta)) ^ (1 / (1 - aalpha))) * (z_avg * nSS);
		output = (K ^ aalpha) * ((z_avg * nSS) ^ (1 - aalpha));
		rK = aalpha * m * (output / K);
		w = (1 - aalpha) * m * (output / (z_avg * nSS));
		
		% Compute labor supply stuff
		cchi = w * (1 - ttau) / (nSS ^ pphi);
		ttransfer = z_cutoff * ttau * w * (nSS ^ (1 / pphi));
		
		% Use this to construct bounds on asset grid with amin = -1 * avg yearly income
		amin = 0;
		%amin = -1 * w * (1 - ttau) * z_avg * nSS;
		amax = 150 * w * (1 - ttau) * z_avg * nSS;

		% Asset grid in positive space (non-linearly spaced grid)
		x = linspace(0,1,I-I_neg)';
		coeff_power = .95; power = 8;
		xx = (1 - coeff_power) * x + coeff_power * (x .^ power);
		xmax = max(xx); xmin = min(xx);
		a_pos = (amax / (xmax - xmin)) * xx + 0;
		clear x xx
		
		% Asset grid in negative space (linearly spaced grid)
		x = linspace(0,1,I_neg+1)';
		coeff_power = 0; coeff_lin = 0; power = 1;
		xx = (1 - coeff_power) * x + coeff_power * (x .^ power);
		xmax = max(xx); xmin = min(xx);
		a_neg = ((0 - amin) / (xmax - xmin)) * xx + amin;
		a_neg(end) = [];
		
		% Put asset grid together
		a = [a_neg;a_pos];
		aa = a * ones(1,J);
		
		% Create grids for finite differences scheme
		daf = ones(I,1);
		dab = ones(I,1);
		daf(1:I-1) = a(2:I) - a(1:I-1);
		dab(2:I) = a(2:I) - a(1:I-1);
		daf(I) = daf(I-1); dab(1) = dab(2);
		daa_f = daf * ones(1,J);
		daa_b = dab * ones(1,J);
		
		% Effective interest rate taking into wedge
		r = ddeath + r_guess .* (aa >= 0) + (r_guess + wedge) .* (aa < 0);
	
		% Initial guess of value function
		v0 = ((w * (1 - ttau) * zz * nSS + r .* aa + ttransfer) .^ (1 - ggamma)) / (rrho_effective * (1 - ggamma));
		disutil = cchi * zz .* ((nSS ^ (1 + pphi)) / (1 + pphi));
		
		% Create initial guess
		if ir>1
			v0 = V_r(:,:,ir-1);
		end
		v = v0;
		
		%%%%
		% Solve for value function given r
		%%%%
		
		for n = 1:maxit
		
			V = v;
			V_n(:,:,n) = V;
			
			% Compute forward difference
			dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) ./ (aa(2:I,:) - aa(1:I-1,:));
			dVf(I,:) = (w * (1 - ttau) * z * nSS + r(I,:) .* amax + ttransfer - disutil(I,:)) .^ (-ggamma);	% will never be used, but impose a <= amax just in case
			
			% Compute backward difference
			dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) ./ (aa(2:I,:) - aa(1:I-1,:));
			dVb(1,:) = (w * (1 - ttau) * z * nSS + r(1,:) .* amin + ttransfer - disutil(1,:)) .^ (-ggamma);	% state constraint boundary condition again
			
			% Consumption and savings with forward difference
			cf = max(dVf,10^(-8)) .^ (-1 / ggamma) + disutil;
			cf(I,:) = min(cf(I,:),w*(1 - ttau)*zz(I,:)*nSS + r(I,:).*aa(I,:) + ttransfer);
			sf = w*(1-ttau)*zz*nSS + r.*aa + ttransfer - cf;
			
			% Consumption and savings with backward difference
			cb = max(dVb,10^(-8)) .^ (-1 / ggamma) + disutil;
			%cb(1,:) = max(cb(1,:),w*(1-ttau)*zz(1,:)*nSS + r(1,:).*aa(1,:)); % Ben commented this out
			sb = w*(1 - ttau)*zz*nSS + r.*aa + ttransfer - cb;
			sb(1,:) = 0;
			
			% Consumption with no drift
			c0 = w*(1 - ttau) * zz*nSS + r.*aa + ttransfer;
			dV0 = (c0 - disutil) .^ (-ggamma);
			
			% dV_Upwind makes choice between forward or backward difference based on sign of drift
			If = sf > 0;	% positive drift --> forward difference
			Ib = (sb < 0) .* (1 - If);	% negative drift --> backward difference
			I0 = 1 - If - Ib; % no drift
			
			% Consumption choice with upwind differences
			c = If .* cf + Ib .* cb + I0 .* c0;
			u = ((c - disutil) .^ (1 - ggamma)) ./ (1 - ggamma);
			dV_Upwind = dVf .* If + dVb .* Ib + dV0 .* I0;
			savingsSS = w*(1-ttau)*zz.*nSS + r.*aa + ttransfer - c;
			
			% Construct matrix for implicit updating scheme
			X = -min(sb,0) ./ daa_b;
			Y = -max(sf,0) ./ daa_f + min(sb,0) ./ daa_b;
			Z = max(sf,0) ./ daa_f;
			
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
			if max(abs(sum(A,2)))>10^(-8)
			   disp('Improper Transition Matrix')
			   break
			end
			
			% Solve for new value function
			B = (1/Delta + rrho_effective)*speye(I*J) - A;
		
			u_stacked = reshape(u,I*J,1);
			V_stacked = reshape(V,I*J,1);
			vec = u_stacked + V_stacked/Delta;
			
			V_stacked = B\vec; % Implicit scheme
			
			V = reshape(V_stacked,I,J);
			Vchange = V - v;
			v = V; 
			
			% Update convergence criterion
			dist(n) = max(max(abs(Vchange)));
			if dist(n)<crit
				%disp('Value Function Converged, Iteration = ')
				%disp(n)
				break
				
			end
			
		end
		
		V_r(:,:,ir) = V;
		
		%%%%
		% Solve for stationary distribution
		%%%%
		
		% Adjust to ensure mass preservation with non-uniform grid
		da_tilde = 0.5*(dab + daf);
		da_tilde(1) = 0.5*daf(1); da_tilde(I) = 0.5*dab(I);
		da_stacked = reshape(da_tilde*ones(1,J),I*J,1);
		grid_diag = spdiags(da_stacked,0,I*J,I*J);
		A_adj = grid_diag*A*grid_diag^(-1); %mass preservation

		% Normalization so pdf integrates to 1 
		gg0 = ones(I*J,1);
		g_sum = gg0'*ones(I*J,1)*da_stacked;
		gg0 = gg0./g_sum;
		
		% Adjust for death
		aux = sparse(I,I);
		aux(:,1)=ddeath*da_tilde/da_tilde(1);
		aux2 = kron(speye(J),aux);
		A_adj = A_adj + aux2 - ddeath*speye(I*J);
		
		% Solve linear system for distribution using implicit method
		AT = A_adj';
		gg = gg0;
		for n = 1:50
			gg_new = (speye(I*J) - AT*Delta)\gg;
			g_sum_t(n) = gg'*da_stacked;
			dist(n) = max(abs(gg_new-gg));
			if dist(n)<10^(-6)
				break
			end
			gg = gg_new;
		end
		
		% Store distribution
		gg_st = gg;
		A_st = (reshape(aa,I*J,1).*gg_st)'*da_stacked;
		C_st = (reshape(c,I*J,1).*gg_st)'*da_stacked;
		c_st = c;
		adot_st = w*(1-ttau)*nSS*zz + r.*aa + ttransfer - c;
		clear gg

		% Marginal distributions
		g = reshape(gg_st,I,J);
		g_a = sum(g,2);
		g_z = da_tilde'*g;
		
		%%%%
		% Update interest rate
		%%%%
		
		S(ir,1) = A_st - K;
        S(ir)
		
		if S(ir)>crit_S
    
			disp('Excess Supply')
			rmax = r_guess;
			r_guess = 0.5*(r_guess+rmin);
			
		elseif S(ir)<-crit_S;
		
			disp('Excess Demand')
			rmin = r_guess;
			r_guess = 0.5*(r_guess+rmax);
			
		elseif abs(S(ir))<crit_S;
		
			display('Steady State Found, Interest rate =')
			disp(r_r(ir))
			break
			
		end
		
	end
	
	%%%%
	% Update discount rate
	%%%%
	%{
	objective(iRrho) = A_st / (w * (1 - ttau) * nSS * z_avg) - 8.85;
	
	if objective(iRrho) > crit_objective
		display('Wealth too high')
		rrhoMin = rrho;
		rrho = 0.5 * (rrho + rrhoMax);
	elseif objective(iRrho) < -crit_objective
		display('Wealth too low')
		rrhoMax = rrho;
		rrho = 0.5 * (rrho + rrhoMin);
	elseif abs(objective(ir)) < crit_objective
		display('Calibration converged')
		display(['Annual discount rate = ' num2str(4*rrho)])
		display(['Capital to output ratio = ' num2str(A_st / output)])
		break
	end
	%}
end
