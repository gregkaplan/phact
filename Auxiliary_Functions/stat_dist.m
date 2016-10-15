function stat_dist = stat_dist(lambda)

eigvl = eig(lambda); [eigvc,~] = eig(lambda);
[~,pos] = min(abs(eigvl));

stat_dist = eigvc(:,pos); stat_dist = stat_dist./sum(stat_dist);

end