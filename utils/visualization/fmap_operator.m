function op = fmap_operator(numEig, L, M)

	[ED, EV] = eigsReal(L, M, numEig);
	EV = massOrthogonal(EV, M);
	op = @compute;

	function [C, nC, nL, nD] = compute(P, Lc, Mc)
		[EDc, EVc] = eigsReal(Lc, Mc, numEig);
		EVc = massOrthogonal(EVc, Mc);
		C = EVc' * Mc * P * EV;
		[nC, nL, nD] = fmap_norms(C, ED, EDc);
	end
end
