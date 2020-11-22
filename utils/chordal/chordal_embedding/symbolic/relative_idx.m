function [relptr,relidx] = relative_idx(colptr, rowidx, snptr, snpar)
%     Compute relative indices of update matrices in frontal matrix of parent.
    
    relptr = ones(length(snptr),1);
    relidx = zeros(colptr(end)-1,1);

    function ret = lfind(a,b)
        i = 1;
        ret = a;
        for k2 = 1:length(a)
            while a(k2) ~= b(i)
                i = i + 1;
            end
            ret(k2) = i;
            i = i + 1;
        end
    end
    for k = 1:length(snpar)
        p = snpar(k);
        relptr(k+1) = relptr(k);
        if p ~= 0
            nk = snptr(k+1) - snptr(k);
            relptr(k+1) = relptr(k+1) + colptr(k+1) - colptr(k) - nk;
            relidx(relptr(k):relptr(k+1)) = lfind(rowidx(colptr(k)+nk:colptr(k+1)), rowidx(colptr(p):colptr(p+1)));
        end
    end
    relidx = relidx(1:relptr(k+1));
end