function [l_grid,w,r_a_grid] = derive_aggregates(KL)

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T

% wage and illiquid asset return

w = (1-alpha) * exp(mean_tfp) * KL^alpha;

r_a = alpha * exp(mean_tfp) * KL^(alpha-1) - delta;
r_a_grid = repmat(r_a,I,J,N);

% unweighted labor input (here trivial)

l = ones(1,N);
l_grid = permute(repmat(l',1,I,J),[2 3 1]);

end