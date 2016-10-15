function x = vertcat(varargin)
% by SeHyoun Ahn, Jan 2016

i = 1;
while (~isa(varargin{i}, 'myAD'))
    i=i+1;
end
y=vertcat(varargin{1:i-1});

x = varargin{i};
l = size(x.derivatives,2);
if (i>1)
    aux = size(y);
    n = numel(y);
    locs = reshape(1:n,aux);
    aux = size(x.values);
    x.values = [y; x.values];
    x.derivatives = [sparse(n, l); x.derivatives];
    locs = [locs; n+reshape(1:prod(aux),aux)];
    n = n+prod(aux);
else
    aux = size(x.values);
    locs = reshape(1:prod(aux),aux);
    n = prod(aux);
end

for j = i+1:nargin
    if isa(varargin{j}, 'myAD')
        x.values = [x.values; varargin{j}.values];
        x.derivatives = [x.derivatives; varargin{j}.derivatives];
        aux = size(varargin{j}.values);
        locs = [locs; n+reshape(1:prod(aux),aux)];
        n = n+prod(aux);
    elseif (~isempty(varargin{j}))
        x.values = [x.values; varargin{j}];
        x.derivatives = [x.derivatives;sparse(numel(varargin{j}),l)];
        aux = size(varargin{j});
        locs = [locs; n+reshape(1:prod(aux),aux)];
        n = n+prod(aux);
    end
end
x.derivatives = x.derivatives(locs(:),:);
