K = N;
a_frac = 1;
b_frac = 1;

b_grid = permute(repmat(b,1,J,K),[1 2 3]);
a_grid = permute(repmat(a,1,I,K),[2 1 3]);
y_grid = permute(repmat(y',1,I,J),[2 3 1]);
r_b_grid = r_b.*(b_grid>=0) + r_b_borr.*(b_grid<0);
T_grid = T * ones(I,J,K);

dbf_grid = NaN(I,J,K);
dbf_grid(1:I-1,:,:) = b_grid(2:I,:,:) - b_grid(1:I-1,:,:); dbf_grid(I,:,:) = dbf_grid(I-1,:,:);
dbb_grid = NaN(I,J,K);
dbb_grid(2:I,:,:) = b_grid(2:I,:,:) - b_grid(1:I-1,:,:); dbb_grid(1,:,:) = dbb_grid(2,:,:);

daf_grid = NaN(I,J,K);
daf_grid(:,1:J-1,:) = a_grid(:,2:J,:) - a_grid(:,1:J-1,:); daf_grid(:,J,:) = daf_grid(:,J-1,:);
dab_grid = NaN(I,J,K);
dab_grid(:,2:J,:) = a_grid(:,2:J,:) - a_grid(:,1:J-1,:); dab_grid(:,1,:) = dab_grid(:,2,:);

db_tilde = 0.5*(dbb_grid(:,1,1) + dbf_grid(:,1,1)); db_tilde(1) = 0.5*dbf_grid(1,1,1); db_tilde(end) = 0.5*dbb_grid(end,1,1);
da_tilde = 0.5*(dab_grid(1,:,1) + daf_grid(1,:,1))'; da_tilde(1) = 0.5*daf_grid(1,1,1); da_tilde(end) = 0.5*dab_grid(1,end,1);
dab_tilde = kron(da_tilde,db_tilde);

D = spdiags(repmat(dab_tilde,K,1),0,I*J*K,I*J*K);
D_inv = spdiags(repmat(ones(I*J,1)./dab_tilde,K,1),0,I*J*K,I*J*K);
dab_tilde_grid = reshape(repmat(dab_tilde,K,1),I,J,K);

gplot = gSS .* dab_tilde_grid;

figure(1)

for k = 1:N

subplot(1,2,k)
set(gca,'FontSize',16)
surf(b(1:round(b_frac*I),1),a(1:round(a_frac*J),1),gplot(1:round(b_frac*I),1:round(a_frac*J),k)')
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
title([ 'Distribution, Type' num2str(k) ])

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 2*pos(4)]);