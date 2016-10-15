X(2:I,:,:) = - bdriftB(2:I,:,:)./dbb_grid(2:I,:,:); X(1,:,:) = zeros(1,J,N);
Y(2:I-1,:,:) = bdriftB(2:I-1,:,:)./dbb_grid(2:I-1,:,:) - bdriftF(2:I-1,:,:)./dbf_grid(2:I-1,:,:); Y(1,:,:) = - bdriftF(1,:,:)./dbf_grid(1,:,:); Y(I,:,:) = bdriftB(I,:,:)./dbb_grid(I,:,:);
Z(1:I-1,:,:) = bdriftF(1:I-1,:,:)./dbf_grid(1:I-1,:,:); Z(I,:,:) = zeros(1,J,N);

centdiag = reshape(Y,I*J,N);
lowdiag = reshape(X,I*J,N);
lowdiag = circshift(lowdiag,-1);
updiag = reshape(Z,I*J,N);
updiag = circshift(updiag,1);

centdiag = reshape(centdiag,I*J*N,1);
updiag   = reshape(updiag,I*J*N,1);
lowdiag  = reshape(lowdiag,I*J*N,1);

bb = spdiags(centdiag,0,I*J*N,I*J*N) + spdiags(updiag,1,I*J*N,I*J*N) + spdiags(lowdiag,-1,I*J*N,I*J*N);

Xu(2:I,:,:) = - budriftB(2:I,:,:)./dbb_grid(2:I,:,:); Xu(1,:,:) = zeros(1,J,N);
Yu(2:I-1,:,:) = budriftB(2:I-1,:,:)./dbb_grid(2:I-1,:,:) - budriftF(2:I-1,:,:)./dbf_grid(2:I-1,:,:); Yu(1,:,:) = - budriftF(1,:,:)./dbf_grid(1,:,:); Yu(I,:,:) = budriftB(I,:,:)./dbb_grid(I,:,:);
Zu(1:I-1,:,:) = budriftF(1:I-1,:,:)./dbf_grid(1:I-1,:,:); Zu(I,:,:) = zeros(1,J,N);

centdiagu = reshape(Yu,I*J,N);
lowdiagu = reshape(Xu,I*J,N);
lowdiagu = circshift(lowdiagu,-1);
updiagu = reshape(Zu,I*J,N);
updiagu = circshift(updiagu,1);

centdiagu = reshape(centdiagu,I*J*N,1);
updiagu   = reshape(updiagu,I*J*N,1);
lowdiagu  = reshape(lowdiagu,I*J*N,1);

bbu = spdiags(centdiagu,0,I*J*N,I*J*N) + spdiags(updiagu,1,I*J*N,I*J*N) + spdiags(lowdiagu,-1,I*J*N,I*J*N);