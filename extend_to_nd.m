function [from_approx, to_approx] = extend_to_nd(from_small,to_small,n_prior,n_post)
% Given a 1-dimensional approximation, approximate n-dimensional function by
%    approximating one of the dimensions using the given 1-dim approximation
%
% by SeHyoun Ahn, June 2016
%
% PARAMETERS:
%    from_small = projection to approximation
%    to_small = projection back from approximation
%    n_prior = number of grid points stacked prior to reduction dimension
%    n_post = number of grid points stacked after reduction dimension
%
% OUTPUTS:
%    from_approx = basis change from small points to large
%    to_approx = basis change from large number of points to small
%
% EXAMPLE:
%    x = linspace(-1,1,100)';
%    y = linspace(-1,1,30);
%    [from_cheb, to_cheb] = cheb_1_dim(x,3);
%    [from_approx, to_approx] = extend_to_nd(from_cheb,to_cheb,1,30);
%    z = bsxfun(@times,exp(x),y.^2);
%    surf(reshape(from_approx*to_approx*z(:),100,30));
%    surf(z);
%
% NOTE: Though it is written as 1 dimension to n dimension, if the approximation
%    is over m-dimensional, it can still be extended over (n-m) dimension
%
% SYNTAX:
% [from_approx, to_approx] = extend_to_nd(from_small,to_small,n_prior,n_post)


from_approx = kron(kron(speye(n_post),from_small),speye(n_prior));
to_approx = kron(kron(speye(n_post),to_small),speye(n_prior));
