function [basis,Q,error] = deflated_block_arnoldi(A,B,m,varargin)
% Block Arnoldi Process with deflation to compute orthnormal basis for
%    Krylov subspace K_m(A,B)
%
% by SeHyoun Ahn, Jan 2017
%
% REFERENCES: Freund, Roland W. "Model reduction methods based on Krylov
%              subspaces." Acta Numerica 12 (2003): 267-319.
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
% [basis,Q,error] = deflated_block_arnoldi(A,B,m,varargin)


if nargin == 5
    basis = varargin{1};
    Q = varargin{2};
else
    [Q,~] = qr(B,0);
    basis = [];
end

for i = 1:(m-1)

    % Implement modified gram-schmidt
    basis = [basis,Q];
    if isa(A,'function_handle')
        aux = A(Q);
    else
        aux = A*Q;
    end
    for j = 1:size(basis,2)
        aux = aux - bsxfun(@times,basis(:,j),(basis(:,j)'*aux));
    end

    % Check for potential deflation
    Q = [];
    for j = 1:size(aux,2)
        weight = sqrt(sum(aux(:,j).^2));
        if weight > sqrt(eps)
            Q = [Q, aux(:,j)/weight];
            for k = j+1:size(aux,2)
                aux(:,k) = aux(:,k) - Q(:,end)*(Q(:,end)'*aux(:,k));
            end
        else
            % Uncomment if you want messages when deflation happens
            % disp('<deflated_block_arnoldi>: Linear dependence, deflating one vector');
        end
    end

    % Second run of modified gram-schmidt to reorthogonalize to reduce
    %    rounding error. Might not be necessary, but it is a safety
    %    measure.
    for j = 1:size(basis,2)
        Q = Q - bsxfun(@times,basis(:,j),(basis(:,j)'*Q));
    end
    Q = bsxfun(@rdivide,Q,sqrt(sum(Q.^2)));
end

if (m == 1)
    basis = Q;
end

if isa(A,'function_handle')
    [~,error] = qr(A(Q) - basis*(basis'*(A(Q))),0);
else
    [~,error] = qr(A*Q - basis*(basis'*(A*Q)),0);
end
