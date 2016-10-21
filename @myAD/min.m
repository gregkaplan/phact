function varargout = min(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if nargin > 2
    if isa(varargin{2},'double')
        direction = varargin{2};
        if direction == 1
            z=x;
            [n,m]=size(x.values);
            [z.values,idx]=min(x.values,varargin{:});
            z.derivatives = x.derivatives(n*(0:(m-1))+idx,:);
            varargout{2} = idx;
        else
            z=x;
            [n,~]=size(x.values);
            [z.values,idx]=min(x.values,varargin{:});
            z.derivatives = x.derivatives((idx-1)*n+(1:n)',:);
            varargout{2} = idx;
        end
    else
        if isa(varargin{1},'myAD')
            idx = x.values>varargin{1}.values;
            z=x;
            z.values = x.values.*(1-idx)+idx.*varargin{1}.values;
            z.derivatives(idx(:),:) = varargin{1}.derivatives(idx(:),:);
            warning('There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.');
        else
            if isempty(varargin{1})
                z=x;
                [n,m]=size(x.values);
                [z.values,idx]=min(x.values,varargin{:});
                z.derivatives = x.derivatives(n*(0:(m-1))+idx,:);
                varargout{2} = idx;
            else
                idx = x.values>varargin{1};
                z=x;
                z.values = x.values.*(1-idx)+idx.*varargin{1};
                z.derivatives(idx(:),:) = 0;
                warning('There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.');
            end
        end
    end
elseif nargin==2
    if isa(varargin{1},'myAD')
        idx = x.values>varargin{1}.values;
        z=x;
        z.values = x.values.*(1-idx)+idx.*varargin{1}.values;
        z.derivatives(idx(:),:) = varargin{1}.derivatives(idx(:),:);
        warning('There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.');
    else
        idx = x.values>varargin{1};
        z=x;
        z.values = x.values.*(1-idx)+idx.*varargin{1};
        z.derivatives(idx(:),:) = 0;
        warning('There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.');
    end
else
    z=x;
    [n,m]=size(x.values);
    [z.values,idx]=min(x.values,varargin{:});
    z.derivatives = x.derivatives(n*(0:(m-1))+idx,:);
    varargout{2} = idx;
end
varargout{1} = z;
