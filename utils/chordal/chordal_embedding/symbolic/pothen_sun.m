
function [snpar, flag] = pothen_sun(par, post, colcount)
%     Find supernodes and supernodal etree.
% 
%     ARGUMENTS
%     par       parent array
% 
%     post      array with post ordering
% 
%     colcount  array with column counts
% 
%     RETURNS
%     snpar     supernodal parent structure 
% 
%     flag      integer vector of length n; if flag[i] < 0, then -flag[i]
%               is the degree of the supernode with repr. vertex i; if
%               flag[i] >= 0, then flag[i] is the repr. vertex to which
%               node i belongs.
    
    n = length(par);
    flag = zeros(n,1);
    snpar = zeros(n,1);
    snodes = n;
%     ch = {}
    for k = 1:n
        ch(k) = {[]};
    end
    for jj = 1:length(post)
        j = post(jj);
        if par(j) > 0
            ch(par(j)) = {[ch{par(j)},j]};
        end
        mdeg = colcount(j) - 1;

        if par(j) ~= 0
            if (mdeg == colcount(par(j))) && (flag(par(j)) == 0)
                % par[j] not assigned to supernode
                snodes = snodes - 1;
                if flag(j) <= 0   % j is a repr. vertex
                    flag(par(j)) = j;
                    flag(j) = flag(j) - 1;
                else             % j is not a repr. vertex
                    flag(par(j)) = flag(j);
                    flag(flag(j)) = flag(flag(j)) - 1;
                end
            end
        else
            if flag(j) <= 0
                snpar(j) = 0;
            else
                snpar(flag(j)) = 0;
            end
        end

        if flag(j) <= 0
            k = j;
        else
            k = flag(j);
        end

        chj = ch{j};
        for ii = 1:length(chj)
            i = chj(ii);
            if flag(i) <= 0
                l = i;
            else
                l = flag(i);
            end
            if  l ~= k
                snpar(l) = k;
            end
        end
    end
    

    repr = find(flag <= 0);
    
    %deg = -flag(repr);

    % renumber etree with number of supernodes
    sn = zeros(n+1,1);
    for k = 1:length(repr)
        sn(repr(k)) = k;
    end
    
    index = snpar(repr);
    index(index==0) = n+1;
    
    snpar = sn(index);
end
