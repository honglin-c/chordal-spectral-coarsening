function [k,left] = lmerge(left, right, offsetl, offsetr, nl, nr)
    
    tmp = ones(nl+nr,1);
    il = 0;
    ir = 0;
    k = 1;
    while ((il < nl) && (ir < nr))
        if (left(offsetl+il) < right(offsetr+ir))
            tmp(k) = left(offsetl+il);
            il = il +  1;
        elseif (left(offsetl+il) > right(offsetr+ir))
          tmp(k) = right(offsetr+ir);
          ir = ir + 1;
        
        else
            tmp(k) = left(offsetl+il);
            il = il + 1;
            ir = ir + 1;
        end
        k = k + 1;
    end
    if (il < nl)
        for i = 0:(nl-il-1)
            tmp(k+i) = left(offsetl+il+i);
        end
        k = k + nl-il;
    end
    if (ir < nr)
        for i = 0:(nr-ir-1)
            tmp(k+i) = right(offsetr+ir+i);
        end
        k = k + nr-ir;
    end
    
    k = k - 1;
    for i = 0:(k-1)
        left(offsetl+i) = tmp(i+1);
    end
    
    
end