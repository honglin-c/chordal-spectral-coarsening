function [colptr, rowidx] = embed(A, colcount, snode, snptr, snpar, snpost)
%     Compute filled pattern.
% 
%        colptr, rowidx = embed(A, colcount, snode, snptr, snpar, snpost)
% 
%     PURPOSE
%     Computes rowindices and column pointer for representative vertices in supernodes.
% 
%     ARGUMENTS
%     A         sparse matrix
% 
%     colcount  vector with column counts
% 
%     snode     vector with supernodes
%  
%     snptr     vector with offsets
% 
%     snpar     vector with supernodal parent indices
% 
%     snpost    vector with supernodal post ordering
% 
%     RETURNS
%     colptr    vector with offsets 
% 
%     rowidx    vector with rowindices 

    Alo = tril_ccs(A);
    cp = Alo.colptr; ri = Alo.rowidx;
    
    N = length(snpar);

    % colptr for compressed cholesky factor
    colptr = ones(N+1,1);
    for k = 1:N
        colptr(k+1) = colptr(k) + colcount(snode(snptr(k)));
    end
    rowidx = zeros(colptr(end)-1,1);
    cnnz = zeros(N,1);

    % compute compressed sparse representation
    for k = 1:N
        p = snptr(k);
        Nk = snptr(k+1)-p;
        nk = cp(snode(p)+1) - cp(snode(p));
        rowidx(colptr(k):colptr(k)+nk-1) = ri(cp(snode(p)):cp(snode(p)+1)-1);
        cnnz(k) = nk;
        for i = 1:(Nk-1)
            nk = cp(snode(p+i)+1)-cp(snode(p+i));
            [cnnz(k),rowidx] = lmerge(rowidx, ri, colptr(k), cp(snode(p+i)), cnnz(k), nk);
        end
    end
    for kk = 1:length(snpost)
        k = snpost(kk);
        p = snptr(k);
        Nk = snptr(k+1)-p;
        if snpar(k) ~= 0
            [cnnz(snpar(k)),rowidx] = lmerge(rowidx,rowidx,colptr(snpar(k)), colptr(k)+Nk,cnnz(snpar(k)), cnnz(k)-Nk);
        end
    end
