function [from_knots, to_knots] = oneDquad_spline(x,knots)
% Creates a quadratic polynomial based spline basis reduction
%
% by SeHyoun Ahn, June 2016
%
% REFERENCES: To be Written
%
% PARAMETERS:
%    x = fine grid points
%    knots = coarser grid points to reduce to
%
% OUTPUTS:
%    from_knots = change of basis from spline basis to finer grid points
%    to_knots = change of basis from finer grid points to spline basis
%
% EXAMPLE:
%     x = linspace(-1,1,100)';
%     knots = linspace(-1,1,3)';
%     [from_knots,to_knots] = oneDquad_spline(x,knots);
%     y = exp(x);
%     plot(x,y,'b-');hold on; plot(x,from_knots*to_knots*y,'r--');
%     legend('Exact Function','Quad Poly Spline Approx');
%
% SYNTAX:
% [from_knots, to_knots] = oneDquad_spline(x,knots)

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
