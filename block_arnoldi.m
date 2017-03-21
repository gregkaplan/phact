function [basis,Q,error] = block_arnoldi(A,B,m,varargin)
% Block Arnoldi Process to Compute orthnormal basis for Krylov subspace
%      K_m(A,B)
%
% by SeHyoun Ahn, Sept 2016
%
% REFERENCES: Jaimoukha, Imad M., and Ebrahim M. Kasenally. "Krylov subspace
%              methods for solving large Lyapunov equations." SIAM Journal 
%              on Numerical Analysis 31.1 (1994): 227-251.
%
% PARAMETERS:
%     A,B,m corresponding to Krylov subspace K_m(A,B)
%     A = (function handle or matrix)
%     B = (matrix)
%     m = (integer)
%
% OUTPUTS:
%    basis = orthonormalized bases vectors for K_m(A,B)
%    Q = Current candidate corresponding to ``A^mB''
%    error = residual for the projection into current subspace
%            needed for error estimating error in <low_rank_lyap.m>
%
% Note: Still uses bsxfun for backward compatibility, but can be
%           updated for explicit expansion for MATLAB 2017a and later
%
% SYNTAX:
% [basis,Q,error] = block_arnoldi(A,B,m,varargin)


if nargin == 5
    basis = varargin{1};
    Q = varargin{2};
else
    [Q,~] = qr(B,0);
    basis = [];
end

for i = 1:(m-1)
    basis = [basis,Q];
    if isa(A,'function_handle')
        [Q,~] = qr(A(Q) - basis*(basis'*(A(Q))),0);
    else
        [Q,~] = qr(A*Q - basis*(basis'*(A*Q)),0);
    end
end
basis = [basis,Q];

% This process loses orthonormality somtimes, so it reorthonormalizes
%    the basis vectors with qr decomposition.
[basis,~] = qr(basis,0);
if isa(A,'function_handle')
    [Q,error] = qr(A(Q) - basis*(basis'*(A(Q))),0);
else
    [Q,error] = qr(A*Q - basis*(basis'*(A*Q)),0);
end
