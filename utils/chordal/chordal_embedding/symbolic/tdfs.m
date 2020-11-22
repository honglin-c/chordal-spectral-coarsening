function [k,head,post,stack] = tdfs(j, k, head, next, post, stack)
%     Depth-first search and postorder of a tree rooted at node j.
    top = 1;
    stack(1) = j;
    
    while (top >= 1)
        p = stack(top);
        i = head(p);
        if i == 0
            top = top - 1;
            post(k) = p;
            k = k + 1;
        else
            head(p) = next(i);
            top = top + 1;
            stack(top) = i;
        end
    end
end