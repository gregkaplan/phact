b_grid = permute(repmat(b,1,J,N),[1 2 3]);
a_grid = permute(repmat(a,1,I,N),[2 1 3]);
y_grid = permute(repmat(y',1,I,J),[2 3 1]);
r_b_grid = r_b.*(b_grid>=0) + r_b_borr.*(b_grid<0);
trans_grid = trans * ones(I,J,N);

dbf_grid = NaN(I,J,N);
dbf_grid(1:I-1,:,:) = b_grid(2:I,:,:) - b_grid(1:I-1,:,:); dbf_grid(I,:,:) = dbf_grid(I-1,:,:);
dbb_grid = NaN(I,J,N);
dbb_grid(2:I,:,:) = b_grid(2:I,:,:) - b_grid(1:I-1,:,:); dbb_grid(1,:,:) = dbb_grid(2,:,:);

daf_grid = NaN(I,J,N);
daf_grid(:,1:J-1,:) = a_grid(:,2:J,:) - a_grid(:,1:J-1,:); daf_grid(:,J,:) = daf_grid(:,J-1,:);
dab_grid = NaN(I,J,N);
dab_grid(:,2:J,:) = a_grid(:,2:J,:) - a_grid(:,1:J-1,:); dab_grid(:,1,:) = dab_grid(:,2,:);

db_tilde = 0.5*(dbb_grid(:,1,1) + dbf_grid(:,1,1)); db_tilde(1) = 0.5*dbf_grid(1,1,1); db_tilde(end) = 0.5*dbb_grid(end,1,1);
da_tilde = 0.5*(dab_grid(1,:,1) + daf_grid(1,:,1))'; da_tilde(1) = 0.5*daf_grid(1,1,1); da_tilde(end) = 0.5*dab_grid(1,end,1);
dab_tilde = kron(da_tilde,db_tilde);
dab_tilde_grid = reshape(repmat(dab_tilde,N,1),I,J,N);
dab_tilde_mat = spdiags(repmat(dab_tilde,N,1),0,I*J*N,I*J*N);

% D = spdiags(repmat(dab_tilde,N,1),0,I*J*N,I*J*N);
% D_inv = spdiags(repmat(ones(I*J,1)./dab_tilde,N,1),0,I*J*N,I*J*N);