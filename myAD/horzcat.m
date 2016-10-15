function x = horzcat(varargin)
% by SeHyoun Ahn, Jan 2016
i = 1;
while (~isa(varargin{i}, 'myAD'))
    i=i+1;
end
y=horzcat(varargin{1:i-1});

x = varargin{i};
if (i>1)
    aux = size(y);
    n = numel(y);
    locs = reshape(1:n,aux);
    x.values = [y, x.values];
    aux = size(x.values);
    x.derivatives = [sparse(n, size(x.derivatives,2)); x.derivatives];
    locs = [locs, n+reshape(1:prod(aux),aux)];
    n = n+prod(aux);
else
    aux = size(x.values);
    locs = reshape(1:prod(aux),aux);
    n = prod(aux);
end

for j = i+1:nargin
    if isa(varargin{j}, 'myAD')
        x.values = [x.values, varargin{j}.values];
        x.derivatives = [x.derivatives; varargin{j}.derivatives];
        aux = size(varargin{j}.values);
        locs = [locs, n+reshape(1:prod(aux),aux)];
        n = n+prod(aux);
    elseif (~isempty(varargin{j}))
        x.values = [x.values, varargin{j}];
        x.derivatives(end+numel(varargin{j}),end) = 0;
        aux = size(varargin{j});
        locs = [locs, n+reshape(1:prod(aux),aux)];
        n = n+prod(aux);
    end
end
x.derivatives = x.derivatives(locs(:),:);