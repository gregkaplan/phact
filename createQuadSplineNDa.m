function [from_spline, to_spline] = createQuadSplineNDa(x,knots,n_prior,n_post,n_rest)
% Creates a quadratic spline basis reduction
%
% by SeHyoun Ahn, June 2016
%
%% PARAMETERS:
%    x = a vector of fine grid points
%    knots = a vector of coarser grid points
%            (does not have to be uniform)
%    n_prior = size of points stacked prior to reduction dimension
%    n_post = size of points stacked after the dimension reduction
%    n_rest = trailing variables that we do not touch
%
%% OUTPUTS:
%    from_knots = change of basis from spline basis to finer grid points
%    to_knots = change of basis from finer grid points to spline basis

n_a = length(x);
n_knots = length(knots);

% Preallocate auxiliary matrices
first_interp_mat = zeros(n_a,n_knots+1);
aux_mat = zeros(n_a,n_knots);

% Find where each points corresponds to on knots
for i = 1:n_a
    loc = sum(knots <= x(i));
    if loc == n_knots
        loc = n_knots -1;
    end
    first_interp_mat(i,loc) = 1-(x(i)-knots(loc))^2/(knots(loc+1)-knots(loc))^2;
    first_interp_mat(i,loc+1) = (x(i)-knots(loc))^2/(knots(loc+1)-knots(loc))^2;
    aux_mat(i,loc) = (x(i)-knots(loc))-(x(i)-knots(loc))^2/(knots(loc+1)-knots(loc));
end

aux_mat2 = spdiags(ones(n_knots,1),0,n_knots,n_knots)...
    +spdiags(ones(n_knots,1),1,n_knots,n_knots);
aux_mat2(end,end) = 0;
aux_mat2(n_knots,1) = 1;
aux_mat3 = spdiags([-2./diff(knots);0],0,n_knots,n_knots+1)...
    +spdiags([2./diff(knots);1],1,n_knots,n_knots+1);

from_knots = sparse(first_interp_mat)+sparse(aux_mat)*(aux_mat2\aux_mat3);
to_knots = (from_knots'*from_knots)\from_knots'*speye(n_a);

from_spline = from_knots;
to_spline = to_knots;

% from_knots = kron(kron(speye(n_post),from_knots),speye(n_prior));
% to_knots = kron(kron(speye(n_post),to_knots),speye(n_prior));
% 
% n_splined = n_prior*(n_knots+1)*n_post;
% n_v = n_prior*n_a*n_post;
% 
% from_spline = spdiags(ones(n_splined+n_rest,1),n_splined-n_v,n_v+n_rest,n_rest+n_splined);
% to_spline = spdiags(ones(n_splined+n_rest,1),n_v-n_splined,n_rest+n_splined,n_v+n_rest);
% from_spline(1:n_v,1:n_splined) = from_knots;
% to_spline(1:n_splined,1:n_v) = to_knots;

end
