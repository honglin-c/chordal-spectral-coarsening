function [colcount, snode, snptr, snpar, snpost] = ...
    amalgamate(colcount, snode, snptr, snpar, snpost, tfill, tsize)
%     Supernodal amalgamation.
% 
%        colcount, snode, snptr, snpar, snpost = ...
%          amalgamate(colcount, snode, snptr, snpar, snpost, merge_function)
%     
%     PURPOSE
%     Iterates over the clique tree in topological order and greedily
%     merges a supernode with its parent if
% 
%        merge_function(|J_{par(k)}|, |J_k|, |N_{par(k)}|, |N_k|)
% 
%     returns True.
% 
%     ARGUMENTS
%     colcount  vector with column counts
% 
%     snode     vector with supernodes
%  
%     snptr     vector with offsets
% 
%     snpar     vector with supernodal parent indices
% 
%     snpost    vector with supernodal post ordering
% 
%     merge_function
%               function 
% 
%     RETURNS
%     colcount  vector with amalgamated column counts
% 
%     snode     vector with amalgamated supernodes
%  
%     snptr     vector with amalgamated offsets
% 
%     snpar     vector with amalgamated supernodal parent indices
% 
%     snpost    vector with amalgamated supernodal post ordering
    N = length(snpost);
    for k = 1:N
        ch(k) = {[]};
    end
%     ch = {}
    for jj = 1:length(snpost)
        j = snpost(jj);
        if snpar(j) > 0
            ch(snpar(j)) = {[ch{snpar(j)};j]};
        end
    end

    
    for k = 1:N
        snlist(k) = {snode(snptr(k):snptr(k+1)-1)};
    end

    Ns = N;
    for kk = 1:length(snpost)
        k = snpost(kk);
        if snpar(k) ~= 0
            snk = snlist{k};
            colk = colcount(snk(1));
            snp = snlist{snpar(k)};
            colp = colcount(snp(1));
            
            nk = length(snk);
            np = length(snp);
            tomerge = merge_size_fill(colp,colk,np,nk, tsize, tfill);
           
            if tomerge
                % merge supernode k and snpar[k]
                
                snlist(snpar(k)) = {sort([snk;snp])};
              
                
                snlist(k) = {[]};
                snp = snlist{snpar(k)};
                colcount(snp(1)) = colp + nk;
                Ns = Ns - 1;
                
                chk = ch{k};
                for cc = 1:length(chk)
                    c = chk(cc);
                    snpar(c) = snpar(k);
                end
                ch(snpar(k)) = {[ch{snpar(k)};chk]};
                
                snpar(k) = 0;
            end
        end
    end

    offset = 0;
    for i = 1:length(snlist)
        s = snlist{i};
        if ~isempty(s)
            L(offset + 1) = i;
            offset = offset + 1;
        end
    end
    
    for i = 1:length(L)
        snl = snlist{L(i)};
        snptr(i+1) = snptr(i) + length(snl);
        snode(snptr(i):snptr(i+1)-1) = snl;
    end

    snpar = snpar(L);
    for i = 1:length(snpar)
        if snpar(i) ~= 0
            snpar(i) = find(L==snpar(i));
        end
    end
    snpost = post_order(snpar);
end