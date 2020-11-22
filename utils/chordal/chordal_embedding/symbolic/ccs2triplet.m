function [I,J] = ccs2triplet(colptr,rowidx)
    I = rowidx;
    J = 0*I;
    J(colptr(1:end-1)) = 1;
    J = cumsum(J);
    
end