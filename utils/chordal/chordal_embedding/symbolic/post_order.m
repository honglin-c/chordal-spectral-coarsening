function p = post_order(parent)
%     Post order a forest.
    n = length(parent);
    k = 1;

    p = ones(n,1);
    head = zeros(n,1);
    next = ones(n,1);
    stack = ones(n,1);

    for j = (n:-1:1)
        if (parent(j) == 0)
            continue
        end
        next(j) = head(parent(j));
        head(parent(j)) = j;
    end

    for j = 1:n
        if ( parent(j) ~= 0)
            continue
        end
        [k,head,p,stack] = tdfs(j, k, head, next, p, stack);
    end
end
  