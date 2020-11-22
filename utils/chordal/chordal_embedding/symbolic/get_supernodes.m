function [snode, snptr, snpar] = get_supernodes(par, post, colcount)
%     Find supernodes.
% 
%     ARGUMENTS
%     par       parent array
% 
%     post      array with post ordering
% 
%     colcount  array with column counts
% 
%     RETURNS
%     snode     array with supernodes; snode[snptr[k]:snptr[k+1]] contains
%               the indices of supernode k
% 
%     snptr     pointer array; snptr[k] is the index of the representative
%               vertex of supernode k in the snode array
% 
%     snpar     supernodal parent structure 
    [snpar, flag] = pothen_sun(par, post, colcount);
%     save('../matlab_common/symbolic/current2.mat','snpar','par','post','colcount')

    n = length(par);
    N = length(snpar);

    snode = ones(n,1);
    snptr = ones(N+1,1);

    for i = 1:n
        slist(i) = {[]};
    end
    for i = 1:n
        f = flag(i);
        if f <= 0
            slist(i) = {[slist{i},i]};
        else
            slist(f) = {[slist{f},i]};
        end
    end

    k = 0; j = 0;
    for i = 1:length(slist)
        sl = slist{i};
        nsl = length(sl);
        if nsl > 0
            snode(k+1:k+nsl) = sl;
            snptr(j+2) = snptr(j+1) + nsl;
            k = k + nsl;
            j = j + 1;
        end
    end
    
end