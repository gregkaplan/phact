function [from_approx, to_approx] = projection_for_subset(from_small,to_small,n_pre,n_post)
% Create transformations when we only transform subset of the points.
%
% by SeHyoun Ahn, June 2016
%
% REFERENCE: To be Written
%
% PARAMETERS:
%    from_small = (nxk) matrix transformation from approximated k points into full grids
%    to_small = (kxn) matrix transformation into approximation with k points
%    n_pre = number of points preceding the transformation
%    n_post = number of points following the transformation
%
% OUTPUTS:
%    from_approx = transformation from smaller approximation to larger grid points 
%    to_approx = transformation from larger grid points to smaller grid points
%
% EXAMPLE:
%    x = linspace(-1,1,100)';
%    y = x.^2;
%    n_trailing = 3;
%    z = [y;randn(n_trailing,1)];
%    n_cheb = 10;
%    [from_cheb,to_cheb] = cheb_1_dim(x,n_cheb);
%    [from_approx,to_approx] = projection_for_subset(from_cheb,to_cheb,0,n_trailing);
%    z_approxed = from_approx*to_approx*z;
%    disp([z(end-n_trailing+1:end),z_approxed(end-n_trailing+1:end)]);
%    plot(x,z(1:end-n_trailing),'b-');hold on;plot(x,z_approxed(1:end-n_trailing),'r--');
%    legend('exact','approximated');
%
% SYNTAX:
% [from_approx, to_approx] = projection_for_subset(from_small,to_small,n_pre,n_post)

[n_full,n_red] = size(from_small);
from_approx = spdiags(ones(n_red+n_pre,1),0,n_full+n_pre,n_pre+n_red);
to_approx = spdiags(ones(n_red+n_pre,1),0,n_pre+n_red,n_full+n_pre);
from_approx(1+n_pre:n_pre+n_full,n_pre+1:n_pre+n_red) = from_small;
to_approx(n_pre+1:n_pre+n_red,n_pre+1:n_pre+n_full) = to_small;
from_approx(end+1:end+n_post,end+1:end+n_post) = speye(n_post);
to_approx(end+1:end+n_post,end+1:end+n_post) = speye(n_post);
