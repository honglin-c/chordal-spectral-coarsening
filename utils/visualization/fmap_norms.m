function [nC, nL, nD] = fmap_norms(C, ED, EDc)
	numEig = size(C, 1);

	% Norm: Frobenius (squared)
	nC = norm(C, 'fro')^2;

	% Norm: laplacian commutativity (squared)
	nL = norm(C * diag(ED) - diag(EDc) * C, 'fro')^2 / nC;

	% Norm: diagonality (squared)
	D = zeros(numEig, 1);
	for k = 1:numEig
		C_k = C(1:k, 1:k);
		D(k) = norm(C_k' * C_k - speye(k), 'fro')^2;
	end
	nD = sum(D .* ((1:numEig).^-1)');
end
