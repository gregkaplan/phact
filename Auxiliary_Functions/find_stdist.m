function [g,gg] = find_stdist(gg,aau,bbu,dab_tilde_mat,maxit_KFE,Delta_KFE,crit_KFE)

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T

A = aau + bbu;
lambda0 = lambda - diag(diag(lambda));  % transition matrix with diagonal killed
lambda0p = lambda0';

gg = dab_tilde_mat * gg;

for n = 1:maxit_KFE
    
gg1 = NaN(I*J,N);

for k = 1:N
    Ak = A(1+(k-1)*(I*J):k*(I*J),1+(k-1)*(I*J):k*(I*J));
    death_inflow = zeros(I,J);
    death_inflow(b==0,1) = sum(gg(1+(k-1)*(I*J):k*(I*J)));
    death_inflow = reshape(death_inflow,I*J,1);
    gk_sum = sum(repmat(lambda0p(k,:),I*J,1) .* reshape(gg,I*J,N),2);
    gg1(:,k) = (speye(I*J) - Delta_KFE * Ak' - Delta_KFE * (lambda(k,k) - deathrate) *speye(I*J))...
        \(gg(1+(k-1)*(I*J):k*(I*J)) + Delta_KFE*gk_sum + Delta_KFE*deathrate*death_inflow);
end

gg1 = reshape(gg1,I*J*N,1);

gg1_sum = sum(gg1);
gg1 = gg1 ./ gg1_sum;

dist = max(abs(gg1-gg));
if mod(n,100) == 0
%     disp(['Iteration ' int2str(n) ', distance is ' num2str(dist)]);
end
if dist<crit_KFE
    disp(['Distribution Converged, Iteration = ' int2str(n)]);
    break
end
gg = gg1;

end

gg = dab_tilde_mat\gg;

g = reshape(gg,I,J,N);

end