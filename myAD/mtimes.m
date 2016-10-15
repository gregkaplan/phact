function x = mtimes(x, y)
% by SeHyoun Ahn, Jan 2016
% Note that in the original package by Martin Fink, this redirected to element-wise multiplication
if isa(x,'myAD')
    [n,m]=size(x.values);
    if isa(y,'myAD')
        if m>1 && size(y,1)==m
            if size(y,2)>1
                z=x;
                tmp=cell(size(y,2),1);
                for j=1:size(y,2)
                    tmp{j,1}=sparse(x.values)*y.derivatives((j-1)*m+(1:m),:) + matdrivXvecval(x.derivatives,y.values(:,j));
                end
                z.derivatives=cell2mat(tmp);
                z.values=x.values*y.values;
                x=z;
            else
                x.derivatives = sparse(x.values)*y.derivatives + matdrivXvecval(x.derivatives,y.values);
                x.values=x.values*y.values;
            end
        elseif max(m,n)==1
            x.derivatives= bsxfun(@times,x.derivatives, y.values(:))+y.derivatives*x.values;
            x.values=x.values*y.values;
        elseif numel(y.values)==1
            x.derivatives=x.derivatives*y.values+bsxfun(@times,x.values(:),y.derivatives);
            x.values=x.values*y.values;
        else
            error('Check that the dimensions match');
        end
    else
        if m>1 && size(y,1)==m
            if size(y,2)>1
                z=x;
                tmp=cell(size(y,2),1);
                for j=1:size(y,2)
                    tmp{j}=matdrivXvecval(x.derivatives,y(:,j));%z.derivatives((j-1)*n+(1:n),:)=matdrivXvecval(x.derivatives,y(:,j));
                end
                z.derivatives=cell2mat(tmp);
                z.values=x.values*y;
                x=z;
            else
                x.derivatives = matdrivXvecval(x.derivatives,y);
                x.values=x.values*y;
            end
        elseif max(m,n)==1
            x.derivatives= bsxfun(@times,x.derivatives, y(:));
            x.values=x.values*y;
        elseif numel(y)==1
            x.derivatives=x.derivatives*y;
            x.values=x.values*y;
        else
            error('Check that the dimensions match');
        end
    end
else
    [n,m]=size(x);
    if m>1 && size(y,1)==m
        if size(y,2)>1
            z=y;
            tmp=cell(size(y,2),1);
            for j=1:size(y,2)
                tmp{j,1}=x*y.derivatives((j-1)*m+(1:m),:);%z.derivatives((j-1)*n+(1:n),:)=sparse(x)*y.derivatives((j-1)*m+(1:m),:);
            end
            z.derivatives=cell2mat(tmp);
            z.values=x*y.values;
            x=z;
        else
            y.derivatives = sparse(x)*y.derivatives;
            y.values=x*y.values;
            x=y;
        end
    elseif max(m,n)==1
        y.derivatives=y.derivatives*x;
        y.values=x*y.values;
        x=y;
    elseif numel(y.values)==1
        y.derivatives=bsxfun(@times,x(:),y.derivatives);
        y.values=x*y.values;
        x=y;
    else
        error('Check that the dimensions match');
    end
end
