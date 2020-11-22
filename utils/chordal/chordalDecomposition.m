function [A_chordal, cliques, vecCliques] = chordalDecomposition(A, tfill, tsize)

    [Sp, Spfill,Sp_unmerged, p, overlap, veclen, cliques, cliques_unmerged] = chordal_embedding(A,tfill,tsize,0);
    
    % compute the chordal extended A named "A_chord" index "p" and its inverse "invp", 
    % such that Ap = A(p,p) and A = Ap(invp,invp)
    [~, invp] = sort(p);
    A_chordal = Spfill(invp,invp);

    % reindexing cliques so that A_chordal(cliques{k}, cliques{k}) is dense
    for k = 1:length(cliques)
        cliques{k} = p(cliques{k});
    end
    
    % compute vecCliques such that
    % z(vecCliques{k}) = vec(A(cliques{k}, cliques{K}))
    vecCliques = cell(length(cliques), 1);
    currentIdx = 1;
    for k = 1:length(cliques)
        ncl = length(cliques{k});
        endIdx = currentIdx+ (ncl*ncl) - 1;
        vecCliques{k} = [currentIdx: endIdx]';
        currentIdx = endIdx+1;
    end
end