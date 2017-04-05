function [to_fine, from_fine] = twoDquad_spline(x_fine,y_fine,x_knots,y_knots)
% Create 2D quadratic polynomial based spline approximation
%
% by SeHyoun Ahn, July 2016
%
% REFERENCE: To be Written
%
% PARAMTERS:
%    x_fine,y_fine = fine grid points
%    x_knots,y_knots = coarser grid to reduce to
%
% OUTPUTS:
%    to_fine = change of basis from spline basis to finer grid points
%    from_fine = change of basis from finer grid points to spline basis
%
% EXAMPLES:
%   x = linspace(-1,1,100)';
%   y = linspace(-1,1,100)';
%
%   % z = e^(-(x^2+y^2));
%   z = exp(-bsxfun(@plus,x.^2,y'.^2));
%   
%   knot_x = linspace(-1,1,4)';
%   knot_y = linspace(-1,1,4)';
%   
%   [to_fine, from_fine] = twoDquad_spline(x,y,knot_x,knot_y);
%   
%   z_approxed = to_fine*from_fine*z(:);
%   
%   surf(z);
%   hold on;
%   disp('press any key');
%   pause();
%   surf(reshape(z_approxed,100,100));
%
% SYNTAX:
% [to_fine, from_fine] = twoDquad_spline(x_fine,y_fine,x_knots,y_knots)

  
n_x_fine = length(x_fine);
n_y_fine = length(y_fine);

n_x_knots = length(x_knots);
n_y_knots = length(y_knots);
h_x = diff(x_knots);
h_y = diff(y_knots);

% Find the Location of Rectangle
x_locs = sum(bsxfun(@ge,x_fine,x_knots'),2);
x_locs = min(x_locs,n_x_knots-1);
x_adjusted = (x_fine - x_knots(x_locs))./h_x(x_locs);

y_locs = sum(bsxfun(@ge,y_fine,y_knots'),2);
y_locs = min(y_locs,n_y_knots-1);
y_adjusted = (y_fine - y_knots(y_locs))./h_y(y_locs);

xx = bsxfun(@times,x_adjusted,ones(1,n_y_fine));
yy = bsxfun(@times,ones(n_x_fine,1),y_adjusted');

xx_locs = bsxfun(@times,x_locs,ones(1,n_y_fine));
yy_locs = bsxfun(@times,ones(n_x_fine,1),y_locs');

% Stacked value of which rectangle the grid point is in
position = xx_locs(:) + (yy_locs(:)-1)*(n_x_knots-1);

% Evaluation of Polynomial Basis
polynomials_stacked = [ones(n_x_fine*n_y_fine,1),  xx(:),  yy(:),  ...
                        xx(:).^2,  xx(:).*yy(:),  yy(:).^2,  ...
                        xx(:).^2.*yy(:),  xx(:).*yy(:).^2];

evaluation = zeros(n_x_fine*n_y_fine,8*(n_x_knots-1)*(n_y_knots-1));

for i=(n_x_knots-1)*(n_y_knots-1):-1:1
    locs = find(position==i);
    evaluation(locs,(i-1)*8+1:8*i) = polynomials_stacked(locs,:);
end

evaluation = sparse(evaluation);

coeff_change = [1, 0, 0, 0, 0, 0, 0, 0;         % f(0,0)
                1, 1, 0, 1, 0, 0, 0, 0;         % f(1,0)
                1, 0, 1, 0, 0, 1, 0, 0;         % f(0,1)
                1, 1, 1, 1, 1, 1, 1, 1;         % f(1,1)
                0, 1, 0, 0, 0, 0, 0, 0;         % dx(0,0)
                0, 0, 1, 0, 0, 0, 0, 0;         % dy(0,0)
                0, 0, 1, 0, 1, 0, 1, 0;         % dy(1,0)
                0, 1, 0, 0, 1, 0, 0, 1];        % dx(0,1)

revert_change = sparse(inv(coeff_change));
aux = kron(speye(n_y_knots-1),spdiags(ones(n_x_knots-1,1),0,n_x_knots-1,n_x_knots));

%% Unpack Values-Derivatives to Proper Location
f00 = [aux,sparse((n_x_knots-1)*(n_y_knots-1),2*n_x_knots*n_y_knots-n_y_knots)]';
f10 = [sparse((n_x_knots-1)*(n_y_knots-1),1),aux,sparse((n_x_knots-1)*(n_y_knots-1),2*n_x_knots*n_y_knots-n_y_knots-1)]';
f01 = [sparse((n_x_knots-1)*(n_y_knots-1),n_x_knots),aux,sparse((n_x_knots-1)*(n_y_knots-1),2*n_x_knots*n_y_knots-n_x_knots-n_y_knots)]';
f11 = [sparse((n_x_knots-1)*(n_y_knots-1),n_x_knots+1),aux,sparse((n_x_knots-1)*(n_y_knots-1),2*n_x_knots*n_y_knots-n_x_knots-n_y_knots-1)]';
dx00 = [sparse((n_x_knots-1)*(n_y_knots-1),n_x_knots*n_y_knots),spdiags(repmat(h_x,n_y_knots-1,1),0,(n_x_knots-1)*(n_y_knots-1),(n_x_knots-1)*(n_y_knots-1)),sparse((n_x_knots-1)*(n_y_knots-1),n_x_knots*n_y_knots-1)]';
dx01 = [sparse((n_x_knots-1)*(n_y_knots-1),n_x_knots*n_y_knots+n_x_knots-1),spdiags(repmat(h_x,n_y_knots-1,1),0,(n_x_knots-1)*(n_y_knots-1),(n_x_knots-1)*(n_y_knots-1)),sparse((n_x_knots-1)*(n_y_knots-1),n_x_knots*n_y_knots-n_x_knots)]';
aux = kron(spdiags(h_y,0,n_y_knots-1,n_y_knots-1),spdiags(ones(n_x_knots-1,1),0,n_x_knots-1,n_x_knots));
dy00 = [sparse((n_x_knots-1)*(n_y_knots-1),2*n_x_knots*n_y_knots-n_y_knots),aux]';
dy10 = [sparse((n_x_knots-1)*(n_y_knots-1),2*n_x_knots*n_y_knots-n_y_knots+1),aux(:,1:end-1)]';

% Warning: Overwriting aux matrix here
aux = [f00(:),f10(:),f01(:),f11(:),dx00(:),dy00(:),dy10(:),dx01(:)]';

% Change Value/Derivatives to Coefficient Basis
aux = revert_change*aux;

%% A nightmare set of transformations. Good luck!
% These are supposed to reorder things to match location of the basis
% expansion
aux = reshape(aux,8*(3*n_x_knots*n_y_knots-n_x_knots-n_y_knots),(n_x_knots-1)*(n_y_knots-1));
permute = reshape(1:8*(3*n_x_knots*n_y_knots-n_x_knots-n_y_knots),8,3*n_x_knots*n_y_knots-n_x_knots-n_y_knots)';
aux = reshape(aux(permute(:),:),3*n_x_knots*n_y_knots-n_x_knots-n_y_knots,8*(n_x_knots-1)*(n_y_knots-1))';

%% Finally Basis Change
to_fine = evaluation*aux;

% Sometimes less than full basis are required up to machine accuarcy, so reduce with SVD
[u,d,~] = svd(full(to_fine),'econ');
singles = find(diag(d)>d(1,1)*eps*10);

% Basis Change Functions
to_fine = u(:,singles)*d(singles,singles);
% Since we took an SVD, so we can just invert by multiplication
from_fine = (d(singles,singles)\u(:,singles)');
