function d = opt_deposits(Va,Vb,a_grid)

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T
  
indx_0     = ((Va./Vb - 1 - chi0) <= 0) .* ((Va./Vb - 1 + chi0) >= 0);
indx_plus  = ((Va./Vb - 1 - chi0) > 0);
indx_minus = ((Va./Vb - 1 + chi0) < 0);

d = 0*indx_0 ...
    + chi1 * (max(Va./Vb - 1 - chi0,0)).^(1/chi2) .* max(a_grid,a_lb) .* indx_plus ...
    + (-chi1) * (max(-(Va./Vb - 1) - chi0,0)).^(1/chi2) .* max(a_grid,a_lb) .* indx_minus;

end