% placeholders

audriftB = NaN(I,J,N);
audriftF = NaN(I,J,N);
budriftB = NaN(I,J,N);
budriftF = NaN(I,J,N);

chi = NaN(I,J,N);
yy = NaN(I,J,N);
zeta = NaN(I,J,N);

X = NaN(I,J,N);
Y = NaN(I,J,N);
Z = NaN(I,J,N);

chiu = NaN(I,J,N);
yyu = NaN(I,J,N);
zetau = NaN(I,J,N);

Xu = NaN(I,J,N);
Yu = NaN(I,J,N);
Zu = NaN(I,J,N);

% compute all drifts
adriftB = min(d,0) + min(a_grid .* (r_a_grid + deathrate*pam) + xi * w * l_grid .* y_grid,0);
adriftF = max(d,0) + max(a_grid .* (r_a_grid + deathrate*pam) + xi * w * l_grid .* y_grid,0);

bdriftB = min(-d - adj_cost_fn(d,a_grid),0) + min(s,0);
bdriftF = max(-d - adj_cost_fn(d,a_grid),0) + max(s,0);

audriftB(1:I-1,:,:) = min(d(1:I-1,:,:) + a_grid(1:I-1,:,:) .* (r_a_grid(1:I-1,:,:) + deathrate*pam) + xi * w * l_grid(1:I-1,:,:) .* y_grid(1:I-1,:,:),0);
audriftB(I,:,:) = min(d(I,:,:) + a_grid(I,:,:) .* (r_a_grid(I,:,:) + deathrate*pam) + xi * w * l_grid(I,:,:) .* y_grid(I,:,:),0);
audriftF(1:I-1,:,:) = max(d(1:I-1,:,:) + a_grid(1:I-1,:,:) .* (r_a_grid(1:I-1,:,:) + deathrate*pam) + xi * w * l_grid(1:I-1,:,:) .* y_grid(1:I-1,:,:),0);
audriftF(I,:,:) = max(d(I,:,:) + a_grid(I,:,:) .* (r_a_grid(I,:,:) + deathrate*pam) + xi * w * l_grid(I,:,:) .* y_grid(I,:,:),0);

budriftB(1:I-1,:,:) = min(s(1:I-1,:,:) - d(1:I-1,:,:) - adj_cost_fn(d(1:I-1,:,:),a_grid(1:I-1,:,:)),0);
budriftB(I,:,:) = min(s(I,:,:) - d(I,:,:) - adj_cost_fn(d(I,:,:),a_grid(I,:,:)),0);
budriftF(1:I-1,:,:) = max(s(1:I-1,:,:) - d(1:I-1,:,:) - adj_cost_fn(d(1:I-1,:,:),a_grid(1:I-1,:,:)),0);
budriftF(I,:,:) = max(s(I,:,:) - d(I,:,:) - adj_cost_fn(d(I,:,:),a_grid(I,:,:)),0);