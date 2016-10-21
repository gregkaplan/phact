function x = subsref(x, S)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

aux=size(x.values);
if ~issparse(x.values)
    locs=reshape(1:prod(aux),aux);
    locs=locs(S.subs{:});
    x.derivatives = x.derivatives(locs(:),:);
    x.values = x.values(S.subs{:});
else
    x.values = x.values(S.subs{:});
    if aux(2)==1
        x.derivatives = x.derivatives(S.subs{1},:);
    elseif aux(1)==1
        x.derivatives = x.derivatives(S.subs{1},:);
    else
        I = aux(1);
        if strcmp(S.subs{1},':')
            S.subs{1} = 1:aux(1);
        end
        if strcmp(S.subs{2},':')
            S.subs{2} = 1:aux(2);
        end

        S.subs{1} = mod(S.subs{1},aux(1));
        func1 = @(i) find(i==S.subs{1});
        func2 = @(i) find(i==S.subs{2});

        n = length(S.subs{1});
        m = length(S.subs{2});
        l = size(x.derivatives,2);
        
        [i,j,v] = find(x.derivatives);
        locs = (ismember(mod(i,I),S.subs{1}) & ismember(ceil(i/I),S.subs{2}));
        i = arrayfun(func1,mod(i(locs),I)) + (arrayfun(func2,ceil(i(locs)/I))-1)*length(S.subs{1});
        j = j(locs);
        v = v(locs);
        x.derivatives = sparse(i,j,v,n*m,l);
    end
end