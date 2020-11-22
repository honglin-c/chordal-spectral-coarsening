function f = objFuncTemplate(var, A, B, L, Mc, invMc)
    f = full(sum(sum(Mc*((A - invMc*var'*L*var*B)).^2))) / 2;
end
