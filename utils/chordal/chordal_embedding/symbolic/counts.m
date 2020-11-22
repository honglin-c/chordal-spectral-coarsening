function colcount = counts(A, parent, post)
  
%     Compute column counts.
    n = A.m;
    colcount = zeros(n,1);
    ancestor = 1:n;
    maxfirst = zeros(n,1);
    prevleaf = zeros(n,1);
    first = zeros(n,1);
    
    for k = 1:n
        j = post(k);
        if first(j) == 0
            colcount(j) = 1;
        else
            colcount(j) = 0;
        end
        while  (j ~= 0) && (first(j) == 0)
            first(j) = k;
            j = parent(j);
        end
    end

    cp = A.colptr; ri = A.rowidx;
    for k = 1:n
        j = post(k);
        if parent(j) ~= 0
            colcount(parent(j)) = colcount(parent(j)) - 1;
        end
        for p = cp(j):(cp(j+1)-1)
            i = ri(p);
            if i <= j
                continue
            end
            
            [q, jleaf, maxfirst, prevleaf,ancestor]  = leaf(i, j, first, maxfirst, prevleaf, ancestor);
            
            if jleaf >= 1
                colcount(j) = colcount(j) + 1;
            end
            if jleaf == 2
                colcount(q) = colcount(q) - 1;
            end
        end
        if parent(j) ~= 0
            ancestor(j) = parent(j);
        end
    end

    for j = 1:n
        if parent(j) ~= 0 
            colcount(parent(j)) = colcount(parent(j)) + colcount(j);
            
        end
    end
end