chi(:,2:J,:) = -adriftB(:,2:J,:)./dab_grid(:,2:J,:); chi(:,1,:) = zeros(I,1,N);
yy(:,2:J-1,:) = adriftB(:,2:J-1,:)./dab_grid(:,2:J-1,:) - adriftF(:,2:J-1,:)./daf_grid(:,2:J-1,:); yy(:,1,:) = - adriftF(:,1,:)./daf_grid(:,1,:); yy(:,J,:) = adriftB(:,J,:)./dab_grid(:,J,:);
zeta(:,1:J-1,:) = adriftF(:,1:J-1,:)./daf_grid(:,1:J-1,:); zeta(:,J,:) = zeros(I,1,N);

centdiag = reshape(yy,I*J,N);
lowdiag = reshape(chi,I*J,N);
lowdiag = circshift(lowdiag,-I);
updiag = reshape(zeta,I*J,N);
updiag = circshift(updiag,I);

centdiag = reshape(centdiag,I*J*N,1);
updiag   = reshape(updiag,I*J*N,1);
lowdiag  = reshape(lowdiag,I*J*N,1);

aa = spdiags(centdiag,0,I*J*N,I*J*N) + spdiags(updiag,I,I*J*N,I*J*N) + spdiags(lowdiag,-I,I*J*N,I*J*N);

chiu(:,2:J,:) = -audriftB(:,2:J,:)./dab_grid(:,2:J,:); chiu(:,1,:) = zeros(I,1,N);
yyu(:,2:J-1,:) = audriftB(:,2:J-1,:)./dab_grid(:,2:J-1,:) - audriftF(:,2:J-1,:)./daf_grid(:,2:J-1,:); yyu(:,1,:) = - audriftF(:,1,:)./daf_grid(:,1,:); yyu(:,J,:) = audriftB(:,J,:)./dab_grid(:,J,:);
zetau(:,1:J-1,:) = audriftF(:,1:J-1,:)./daf_grid(:,1:J-1,:); zetau(:,J,:) = zeros(I,1,N);

centdiagu = reshape(yyu,I*J,N);
lowdiagu = reshape(chiu,I*J,N);
lowdiagu = circshift(lowdiagu,-I);
updiagu = reshape(zetau,I*J,N);
updiagu = circshift(updiagu,I);

centdiagu = reshape(centdiagu,I*J*N,1);
updiagu   = reshape(updiagu,I*J*N,1);
lowdiagu  = reshape(lowdiagu,I*J*N,1);

aau = spdiags(centdiagu,0,I*J*N,I*J*N) + spdiags(updiagu,I,I*J*N,I*J*N) + spdiags(lowdiagu,-I,I*J*N,I*J*N);