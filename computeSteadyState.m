function [raSS,wSS,KSS,uSS,cSS,dSS,VSS,gSS,dab] = computeSteadyState()

%% PREPARATIONS

% global variables

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T

% convergence criteria, size of steps

maxit_HJB = 100; % HJB iterations
crit_HJB = 10^(-5);
Delta = 10^6; % large steps work even though method is partially explicit

maxit_KL = 1000; % loops over capital-labor ratio
crit_KL = 10^(-4);
relax_KL = 0.05;

maxit_KFE = 1000; % KFE iterations
crit_KFE = 10^(-6);
Delta_KFE = 1000; % this is for the implicit method

maxit_HIS = 10; % Howard improvement steps
start_HIS = 2;
crit_HIS = 10^(-5);

% useful large matrices (various large grids, forward and backward steps,
% adjustment matrix for KFE analysis)

set_grids

%% INITIAL GUESSES

% derive all relevant equilibrium aggregates

[l_0_grid,w_0,r_a0_grid] = derive_aggregates(KL_0);

% initial consumption and derived utility guess

c_0 = (1-xi) * w_0 * y_grid .* l_0_grid + (r_a0_grid(1,1,1) + deathrate*pam) .* a_grid ...
    + (r_b_borr + deathrate*pam) .* b_grid + trans_grid - tau_I * w_0 * y_grid .* l_0_grid;
V_0  = 1/(rho + deathrate) * u_fn(c_0,gamma);

% initial distribution guess

gg0 = zeros(I,J,N);
gg0(b==0,1,:) = y_dist;

gg0 = gg0/sum(sum(sum(gg0)));
gg0 = gg0./dab_tilde_grid; % ensure integration to 1

gg0 = reshape(gg0,I*J*N,1);
gg = gg0;

% initialize KL for loop

KL = KL_0;

%% MAIN LOOP

for j = 1:maxit_KL
    
% derive all relevant equilibrium aggregates
    
[l_grid,w,r_a_grid] = derive_aggregates(KL);

% current value function

Vn  = V_0;

%% solve HJB

for n = 1:maxit_HJB
    % derivatives w.r.t. a
    % preparations
    VaF = zeros(I,J,N);
    VaB = zeros(I,J,N);
    Vamin = 0;
    
    % forward difference
    VaF(:,1:J-1,:) = (Vn(:,2:J,:)-Vn(:,1:J-1,:))./daf_grid(:,1:J-1,:);
    VaF(:,1:J-1,:) = max(VaF(:,1:J-1,:),Vamin);
    
    % backward difference
    VaB(:,2:J,:) = (Vn(:,2:J,:)-Vn(:,1:J-1,:))./dab_grid(:,2:J,:);
    VaB(:,2:J,:) = max(VaB(:,2:J,:),Vamin);
    
    % derivatives w.r.t. b
    % preparations
    VbF = zeros(I,J,N);
    VbB = zeros(I,J,N);
    Vbmin = 1e-8;
    
    % forward difference
    VbF(1:I-1,:,:) = (Vn(2:I,:,:)-Vn(1:I-1,:,:))./dbf_grid(1:I-1,:,:);
    VbF(1:I-1,:,:) = max(VbF(1:I-1,:,:),Vbmin);
    
    % backward difference
    VbB(2:I,:,:) = (Vn(2:I,:,:)-Vn(1:I-1,:,:))./dbb_grid(2:I,:,:);
    VbB(2:I,:,:) = max(VbB(2:I,:,:),Vbmin);
    
    % consumption decision  
    % step 0: preparations
    cF = NaN(I,J,N);
    sF = NaN(I,J,N);
    HcF = NaN(I,J,N);
    
    cB = NaN(I,J,N);
    sB = NaN(I,J,N);
    HcB = NaN(I,J,N);
    
    c0 = NaN(I,J,N);
    s0 = NaN(I,J,N);
    Hc0 = NaN(I,J,N);
    
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
    dFB = NaN(I,J,N);
    HdFB = NaN(I,J,N);
    
    dBF = NaN(I,J,N);
    HdBF = NaN(I,J,N);
    
    dBB = NaN(I,J,N);
    HdBB = NaN(I,J,N);
    
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
    transition_prepare    
    
    % transition matrix for a
    transition_a
    
    % transition matrix for b
    transition_b
   
    % transition matrix for productivity state
    transition_c
    
    % overall transition matrix
    
    A = aa + bb;

    % solve for next iteration

    Vn1 = NaN(I,J,N);
    Bik_all = cell(N,1);

    for k = 1:N
        Ak = A(1+(k-1)*(I*J):k*(I*J),1+(k-1)*(I*J):k*(I*J));
        Bk = (1 + Delta*(rho + deathrate) - Delta*lambda(k,k))*speye(I*J) - Delta*Ak;
        Bik_all{k} = inverse(Bk);
        uk_stacked = reshape(u(:,:,k),I*J,1);
        Vk_stacked = reshape(Vn(:,:,k),I*J,1);
        indx_k = ~ismember(1:N,k);
        Vkp_stacked = sum(repmat(lambda(k,indx_k),I*J,1) .* reshape(Vn(:,:,indx_k),I*J,N-1),2);
        qk = Delta*uk_stacked + Vk_stacked + Delta*Vkp_stacked;
        Vn1k_stacked = Bik_all{k}*qk;
        Vn1(:,:,k) = reshape(Vn1k_stacked,I,J,1);
    end
    
    % Howard improvement step
    
    if n >= start_HIS
    for jj = 1:maxit_HIS
    Vn2 = NaN(I,J,N);
    for k = 1:N
        uk_stacked = reshape(u(:,:,k),I*J,1);
        Vk_stacked = reshape(Vn1(:,:,k),I*J,1);
        indx_k = ~ismember(1:N,k);
        Vkp_stacked = sum(repmat(lambda(k,indx_k),I*J,1) .* reshape(Vn1(:,:,indx_k),I*J,N-1),2);
        qk = Delta*uk_stacked + Vk_stacked + Delta*Vkp_stacked;
        Vn2k_stacked = Bik_all{k}*qk;
        Vn2(:,:,k) = reshape(Vn2k_stacked,I,J,1);
    end
    VHIS_delta = Vn2 - Vn1;
    Vn1 = Vn2;
    dist = max(max(max(abs(VHIS_delta))));
    if dist<crit_HIS
        break
    end
    end
    end
    
    % check for convergence
    V_delta = Vn1 - Vn;
    Vn = Vn1;
    dist = max(max(max(abs(V_delta))));
%     disp(['Iteration ' int2str(n) ', distance is ' num2str(dist)]);
        
    if dist<crit_HJB
        disp(['Value Function Converged, Iteration = ' int2str(n)]);
        break
    end
    
end

V_0 = Vn;

%% stationary distribution

[g,gg] = find_stdist_new(gg,aau,bbu,dab_tilde_mat,maxit_KFE,Delta_KFE,crit_KFE);

%% update guess

% compute implied capital supply

K_sup = sum(sum(sum(g.*a_grid.*dab_tilde_grid)));

% compute implied labor supply

L_sup = sum(sum(sum(y_grid .* l_grid .* g .* dab_tilde_grid)));

% compute new capital-labor ratio

KL_sup = K_sup/L_sup;

% check for convergence and update if necessary

gap = KL - KL_sup;
disp(['The current gap is ' num2str(gap)]);

if abs(gap)>crit_KL
    KL = relax_KL*KL_sup + (1-relax_KL)*KL;
else
    disp(['I have found the steady state, Iteration = ' int2str(j)]);
    break
end

end

%% COLLECT FINAL OUTPUTS

KLSS = KL;
raSS = alpha * exp(mean_tfp) * KLSS^(alpha-1) - delta;
wSS  = (1-alpha) * exp(mean_tfp) * KLSS^alpha;
KSS  = KLSS * L_sup;
uSS  = u;
cSS  = c;
dSS  = d;
VSS  = Vn;
gSS  = g;
dab  = dab_tilde_grid;

end