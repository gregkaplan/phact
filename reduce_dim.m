function [dim_red,inv_dim_red] = reduce_dim(V_keep,T,n_v)
% Reduce Dimensionality of State Space
% Written by SeHyoun Ahn, May 2016

n_total = n_v + size(T,2);
n_keep = 0;
V_g = V_keep;
l = size(V_g,2);
for iRun = 1:100
    [~,d_0,v_0] = svds(V_g'*(T/norm(T)),l);
    % The following computes the direction that matter for forecasting
    [V_g,tmp1,~] = svd([V_keep,V_g,bsxfun(@times,v_0,diag(d_0)')],'econ');
    n_old=n_keep;
    %n_keep = sum(diag(tmp1)>10*eps);   % Check the new dimension
	n_keep = sum(diag(tmp1)>sqrt(eps)); 
    V_g = V_g(:,1:n_keep);             % Keep directions with high singular values
    if n_old==n_keep                   % Stop iteration if no new vectors are added.
        n_g=n_keep;
        break;
    end
    tmp1=diag(tmp1);
    V_g=bsxfun(@times,V_g,tmp1(1:n_keep)');
end
n_g = n_keep;

dim_red = sparse(n_v+n_g,n_total);
dim_red(1:n_v,1:n_v)=speye(n_v);
dim_red(n_v+1:end,n_v+1:n_total)=V_g';
inv_dim_red = sparse(n_total,n_v+n_g);
inv_dim_red(1:n_v,1:n_v)=speye(n_v);
inv_dim_red(n_v+1:n_total,n_v+1:n_g+n_v)=V_g;

end
