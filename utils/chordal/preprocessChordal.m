function [X_chordal, cliques, vecCliques, P] = preprocessChordal(X, tfill, tsize)
    % X_chordal: chordal extended X
    % cliques: clique index Zk = A_chordal(clique{k},clique{k})
    % vecCliques: vecclique index Z(vecCliques{k}) = A_chordal(clique{k},clique{k})
    % P: a set of index selection matrix

    if(~exist('tfill','var'))
        tfill = 50;
    end
    if(~exist('tsize','var'))
        tsize = 50;
    end
    
    % check inupt
    assert( all(vec(X-X') == 0)) % assert symmetric
    X = X~=0; % assert binary
    nV = size(X,1);

    % chordal decomposition
    [X_chordal, cliques, vecCliques] = chordalDecomposition(X, tfill, tsize);

    % Compute nonzero selecting matrices such that
    % We define x = P_X2x * X(:) and X(:) = P_x2X * x;
    % where x is the nonzeros of the lower triangular part
    [P.X2x, P.x2X] = selectionMatrixSymmetric(X);
    assert(size(P.X2x,1) == nnz(tril(X)));
    assert(size(P.x2X,2) == nnz(tril(X)));
    
    [P.X2xc, P.xc2X] = selectionMatrixSymmetric(X_chordal);
    assert(size(P.X2xc,1) == nnz(tril(X_chordal)));
    assert(size(P.xc2X,2) == nnz(tril(X_chordal)));
    
    % Compute selection matrices
    % 1. Zk(:) = P_zk2Zk{k} * zk, where zk is the lower triangular part of symmetric Zk
    % 2. zk = P_Zk2zk{k} * Zk(:)
    % 3. X(:) = \sum P_Zk2X{k} * Zk(:)
    P.zk2Zk = cell(length(cliques),1);
    P.Zk2zk = cell(length(cliques),1);
    P.Zk2X = cell(length(cliques),1);
    nz = 0;
    nZ = 0;
    for k = 1:length(cliques)
        if mod(k, 100) == 0
            fprintf('build selection matrices %d / %d\n', k, length(cliques))
        end
        ncl = length(cliques{k});
        Zksp = ones(ncl,ncl);
        [P.Zk2zk{k}, P.zk2Zk{k}] = selectionMatrixSymmetric(Zksp);

        Pk = sparse([1:ncl], cliques{k}, ones(ncl,1), ncl, nV); % X = sum_k Pk'*Zk*Pk
        P.Zk2X{k} = kron(Pk',Pk');

        nz = nz + size(P.Zk2zk{k},1);
        nZ = nZ + size(P.Zk2zk{k},2);
    end
    
    P.z2Z = blkdiag(P.zk2Zk{:});
    P.Z2z = blkdiag(P.Zk2zk{:});
%     assert( all(vec(P.Z2z * P.z2Z - speye(size(P.Z2z,1))) == 0))
    
    % selection matrix between x and xc
    P.x2xc = P.X2xc * P.x2X;
    P.xc2x = P.X2x * P.xc2X;
%     assert( all(vec(P.xc2x * P.x2xc - speye(size(P.xc2x,1))) == 0))
        
    %% other selection matrices
    P.Z2X = horzcat(P.Zk2X{:});
    
    % precomputation
    P.z2X = P.Z2X * P.z2Z;
    P.z2xc = P.X2xc * P.z2X;
    
%     % P.X2Z is ill-defined, people often use simple average
    overlap = sum(P.Z2X,2);
    P.X2Z_avg = P.Z2X';
    [i,j] = find(P.X2Z_avg);
%     idx = find(P.X2Z_avg);
    idx = i + (j-1) * size(P.X2Z_avg,1);
    P.X2Z_avg(idx) = P.X2Z_avg(idx) ./ overlap(j);
    P.X2z_avg = P.Z2z * P.X2Z_avg;
    assert( all(abs(P.X2Z_avg * P.Z2X * ones(size(P.Z2X,2),1) - 1) < 1e-7))
    assert( all(abs(P.X2z_avg * P.z2X * ones(size(P.z2X,2),1) - 1) < 1e-7))

    maxCliqueSize = 0;
    for k = 1:length(cliques)
        if maxCliqueSize < length(cliques{k})
            maxCliqueSize = length(cliques{k});
        end
    end
     
    fprintf("==================================\n")
    fprintf("Chordal decomoposition statistics:\n\n")
    fprintf("number of cliques: %d\n", length(cliques))
    fprintf("maximum clique size: %d\n", maxCliqueSize)
    fprintf("number of variables in X: %d\n", size(P.xc2x,1))
    fprintf("number of variables in Z: %d\n", size(P.z2xc,2))
    fprintf("==================================\n")
    
