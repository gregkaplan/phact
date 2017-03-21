function [from_cheb, to_cheb] = cheb_1_dim(x,n_cheby)
% Create Chebychev polynomial approximation for a function
%
% by SeHyoun Ahn, June 2016
%
% REFERENCES: To be Written
%
% PARAMETERS:
%    x = (nx1 vector) fine grid points (normalized to be between -1 and 1)
%    n_cheby = degree of Chebyshev polynomials to use
%
% OUTPUTS:
%    from_cheb = basis change from polynomial to x-grid values
%    to_cheb = grid value to Chebyshev polynomial values
%
% EXAMPLE:
%     x = linspace(-1,1,100)';
%     [from_cheb,to_cheb] = cheb_1_dim(x,10);
%     y = exp(x);
%     plot(x,y,'b-');hold on; plot(x,from_cheb*to_cheb*y,'r--');
%     legend('Exact Function','Chebyshev Polynomial Approximation');
%
% SYNTAX:
% [from_cheb, to_cheb] = cheb_1_dim(x,n_cheby)


if size(x,2)~=1
    error('Did you mean transpose of x instead of x?')
end

from_cheb=ones(length(x),n_cheby);
from_cheb(:,2)=x;
for i=3:n_cheby
    from_cheb(:,i)=2*x.*from_cheb(:,i-1)-from_cheb(:,i-2);
end
to_cheb = (from_cheb'*from_cheb)\from_cheb'*eye(length(x));
