function L_proj = projectPSD(L,M,epsilon)

    if(~exist('epsilon','var'))
        epsilon = 0;
    end
    
    if max(max(abs(M - speye(size(M,1))))) < 1e-10 
        % mass matrix is identity
        [V,D] = eig(full(L));
        eVal = diag(D);
        eVal(eVal < epsilon) = epsilon;
        D_recon = diag(eVal);
        L_proj = V*D_recon*V';
    else
        % generalized eigenvalue problem 
        [U,D,W] = eig(full(L),full(M));
        assert( abs((U(:,2)' * M * U(:,2))-1) < 1e-7)
        eVal = diag(D);
        eVal(eVal < epsilon) = epsilon;
        D_recon = diag(eVal);
        
        % W'*L = D*W'*M, W'*M*W = I
        % => (W')^(-1) = M*W
        % => L = M*W*D*W'*M
        L_proj = M*W*D_recon*W'*M;
    end
    
    assert(max(max( abs(L_proj - L_proj') )) < 1e-5 )
    L_proj = (L_proj+L_proj') / 2;
end