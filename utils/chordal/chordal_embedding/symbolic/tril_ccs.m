function B = tril_ccs(A)
    colptr = A.colptr;
    I = A.rowidx;
    val = A.val;
    m = A.m;
    n = A.n;
    
    dJ = colptr(2:end-1);
    J = I*0;
    J(dJ) = 1;
    J = cumsum(J)+1;
    
    B = sparse(tril(sparse(I,J,val,m,n)));
    
    [rowidx,colptr,val,n,m] = sparse2ccs(B);
    clear B
    B.rowidx = rowidx;
    B.colptr = colptr;
    B.val = val;
    B.n = n;
    B.m = m;
    
end