function select = select_vec(bstart,bend,astart,aend,ystart,yend)

global rho gamma deathrate chi0 chi1 chi2 a_lb alpha delta rho_tfp sigma_tfp mean_tfp ...
    r_b borrwedge r_b_borr trans tau_I xi pam y y_dist lambda a b I J N KL_0 T ...
    nVars nEErrors rSS wSS KSS uSS cSS dSS VSS gSS dab varsSS

select = zeros(I,J,N);
select(bstart:bend,astart:aend,ystart:yend) = ones(bend-bstart+1,aend-astart+1,yend-ystart+1);
select = reshape(select,I*J*N,1);
select = find(select);

end