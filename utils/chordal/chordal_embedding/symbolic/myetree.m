function parent = myetree(A)
% %     Compute elimination tree from upper triangle of A.
%     assert isinstance(A,spmatrix), "A must be a sparse matrix"
%     assert A.size[0] == A.size[1], "A must be a square matrix"
if A.n ~= A.m
    error('Matrix must be square')
end
    n = A.m;
    cp = A.colptr;
    ri = A.rowidx;
    parent = ones(n,1);
    w = ones(n,1);

    for k = 1:n
        parent(k) = 0;
        w(k) = 0;
        for p = (cp(k)):(cp(k+1)-1)
            i = ri(p);
            while (( i ~= 0) && (i < k))
                inext = w(i);
                w(i) = k;
                if inext == 0 
                    parent(i) = k;
                end
                i = inext;
            end
        end
    end
    
end