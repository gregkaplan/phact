function adj_cost = adj_cost_fn(d,a_grid)

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T

d_scaled = d./max(a_grid,a_lb);
adj_cost = max(a_grid,a_lb) .* (chi0 * abs(d_scaled) + 1/(1+chi2) * (abs(d_scaled).^(1+chi2) * chi1^(-chi2)));

end