function tomerge = merge_size_fill(ccp,cck,np,nk, tsize, tfill)
%     Simple heuristic for supernodal amalgamation (clique
%     merging).
%
%     Returns a function that returns `True` if either (i) supernode k and
%     supernode par(k) are both of size at most `tsize`, or (ii),
%     merging supernodes par[k] and k induces at most `tfill` nonzeros
%     in the lower triangular part of the sparsity pattern.
%
%      PURPOSE
%         Returns true if either (i) supernode k and supernode par(k)
%         are both of size at most %i, or (ii), merging supernodes
%         par(k) and k induces at most %i edges.
%
%         ARGUMENTS
%         Jp        integer; size of parent clique
%
%         Jk        integer; size of child clique
%
%         Np        integer; size of parent supernode
%
%         Nk        integer; size of child supernode
%
%         RETURNS
%         d         boolean

%
%     :param tsize:   nonnegative integer; threshold for merging based on supernode sizes
%     :param tfill:   nonnegative integer; threshold for merging based on induced fill
if tsize < 0
    error('tsize must be nonnegative')
end
if tfill < 0
    error('tfill must be nonnegative')
end


%         Supernodal merge heuristic.
%
%            d = fmerge(Jp, Jk, Np, Nk)
%
%
fill = (ccp - (cck - nk))*nk;
if (fill <= tfill) || ((nk <= tsize) && (np <= tsize))
    tomerge =  true;
else
    tomerge =  false;
end
