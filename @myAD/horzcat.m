function x = horzcat(varargin)
% by SeHyoun Ahn, Jan 2016
i = 1;
while (~isa(varargin{i}, 'myAD'))
    i=i+1;
end
y=horzcat(varargin{1:i-1});

x = varargin{i};
if (i>1)
    n = numel(y);
    x.values = [y, x.values];
    x.derivatives = [sparse(n, size(x.derivatives,2)); x.derivatives];
end

for j = i+1:nargin
    if isa(varargin{j}, 'myAD')
        x.values = [x.values, varargin{j}.values];
        x.derivatives = [x.derivatives; varargin{j}.derivatives];
    elseif (~isempty(varargin{j}))
        x.values = [x.values, varargin{j}];
        x.derivatives(end+numel(varargin{j}),end) = 0;
    end
end
