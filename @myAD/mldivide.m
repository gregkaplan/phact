function y=mldivide(x,y)
% by SeHyoun Ahn, Jan 2016

if isa(x,'myAD')
    [n,m]=size(x.values);
    if isa(y,'myAD')
        if m>1 && size(y,1)==m
            if size(y,2)>1
                z=myAD(x.values\y.values,sparse(n*size(y,2),size(y.derivatives,2)));
                for j=1:size(y,2)
                    z.derivatives((j-1)*n+(1:n),:)=sparse(x.values)\(y.derivatives((j-1)*m+(1:m),:) - matdrivXvecval(x.derivatives,z.values(:,j)));
                end
                y=z;
            else
                y.values = x.values\y.values;
                y.derivatives = sparse(x.values)\(y.derivatives - matdrivXvecval(x.derivatives,y.values));
            end
        elseif max(m,n)==1
            y.derivatives = y.derivatives/x.values - bsxfun(@times,y.values(:)/x.values^2,x.derivatives);
            y.values=y.values/x.values;
        else
            error('Check that the dimensions match');
        end
    else
        if m>1 && size(y,1)==m
            if size(y,2)>1
                z=myAD(x.values\y,sparse(n*size(y,2),size(x.derivatives,2)));
                for j=1:size(y,2)
                    z.derivatives((j-1)*n+(1:n),:)=-sparse(x.values)\(matdrivXvecval(x.derivatives,z.values(:,j)));
                end
                y=z;
            else
                z=myAD(x.values\y,sparse(n*size(y,2),size(x.derivatives,2)));
                z.derivatives = -sparse(x.values)\matdrivXvecval(x.derivatives,z.values);
                y=z;
            end
        elseif max(m,n)==1
            x.derivatives = - bsxfun(@times,y(:)/x.values^2,x.derivatives);
            x.values=y/x.values;
            y=x;
        else
            error('Check that the dimensions match');
        end
    end
else
    [n,m]=size(x);
    if m>1 && size(y,1)==m
        if size(y,2)>1
            z=myAD(x\y.values,sparse(n*size(y,2),size(y.derivatives,2)));
            for j=1:size(y,2)
                z.derivatives((j-1)*n+(1:n),:)=sparse(x)\y.derivatives((j-1)*m+(1:m),:);
            end
            y=z;
        else
            y.values = x\y.values;
            y.derivatives = sparse(x)\y.derivatives;
        end
    elseif max(m,n)==1
        y.derivatives = y.derivatives/x;
        y.values=y.values/x;
    else
        error('Check that the dimensions match');
    end
end
