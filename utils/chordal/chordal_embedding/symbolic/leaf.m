function [q, jleaf, maxfirst, prevleaf,ancestor] = leaf(i, j, first, maxfirst, prevleaf, ancestor)
%     Determine if j is leaf of i'th row subtree.
    jleaf = 0;
    if (i<=j) || (first(j) <= maxfirst(i))
        q = 0;
        return;
    end
    maxfirst(i) = first(j);
    jprev = prevleaf(i);
    prevleaf(i) = j;
    if jprev == 0
        jleaf = 1;
    else
        jleaf = 2;
    end
    if jleaf == 1
        q = i;
        return
    end
    q = jprev;
    while q ~= ancestor(q)
        q = ancestor(q);
    end
    s = jprev;
    while s ~= q
        sparent = ancestor(s);
        ancestor(s) = q;
        s = sparent;
    end
end