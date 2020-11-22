function [P_A2a, P_a2A] = selectionMatrixSymmetric(A)
    % selectionMatrixSymmetric compute selection matrices given sparsity
    % pattern
    %
    % Input:
    %   A is the sparsity pattern of a symmetric matrix 
    % 
    % Output:
    %   P_a2A: selection matrix
    %   P_A2a: selection matrix
    %
    % Note: 
    % - a: a vector non-zeros in the lower triangular part of A
    % - A(:) = P_a2A * a
    % - a = P_A2a * A(:)
    
    n = size(A,1);
    A = A~=0; 
    
    % get P_A2a
    idx = find(tril(A));
    P_A2a = sparse([1:length(idx)], idx, 1, length(idx), n*n);

    % get P_a2A
    tmp = sparse(n,n);
    tmp(idx) = [1:length(idx)];
    tmp = tmp + tmp';
    tmp = tmp - diag(diag(tmp))/2;
    idx_tmp = find(tmp);
    P_a2A = sparse(idx_tmp, tmp(idx_tmp), 1, length(A(:)), size(P_A2a,1));
    
    % check output A(:) == P_a2A * P_A2a * A(:)
    [i,j] = find(A);
    randA = sparse(i,j,rand(length(i),1), n, n);
    randA = (randA + randA');
    assert( max(abs(randA(:) - P_a2A*P_A2a*randA(:))) < 1e-7 )
    na = min(size(P_a2A));
%     assert( all(vec(P_A2a * P_a2A - speye(na)) == 0))
end