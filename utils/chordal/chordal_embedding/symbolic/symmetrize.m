function A = symmetrize(A)
    if A.m ~= A.n
        error('Matrix  must be square')
    end
    [I,J] = ccs2triplet(A.cp,A.ri);
    A2 = sparse(I,J,A.val,A.m,A.n);
    A2= tril(A2) + tril(A2,-1)';
    [rowidx,colptr,val,n,m] = sparse2ccs(A2);
    A.rowidx = rowidx;
    A.colptr = colptr;
    A.val = val;
end