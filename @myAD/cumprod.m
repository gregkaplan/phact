function x = cumprod(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if nargin>1 && varargin{1}==2
    x=ctranspose(cumprod(ctranspose(x)));
else
    [n,m] = size(x.values);
    for i = 2:n
        for j = 1:m
            x.derivatives((j-1)*m+i,:) = x.values(i-1,j).*x.derivatives((j-1)*m+i,:) + x.values(i,j).*x.derivatives((j-1)*m+i-1,:);
        end
        x.values(i,:) = x.values(i-1,:).*x.values(i,:);
    end
end
