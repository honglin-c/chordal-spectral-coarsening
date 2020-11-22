function [rowidx,colptr,val,n,m] = sparse2ccs(A)
    [m,n] = size(A);
    
    
    [rowidx,J,val] = find(sparse(A));
    if length(find(rowidx==J)) < n
        if m == n
            A = A + speye(n);
            [rowidx,J,val] = find(sparse(A));
            val(rowidx==J) = val(rowidx==J) - 1;
        end
    end
    dJ = diff(J);
    [dJI,~,~] = find(sparse(dJ));
    colptr = [1; dJI+1; length(val)+1];
end