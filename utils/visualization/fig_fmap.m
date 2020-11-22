function fig_fmap(fig_id, label, C, nL, nD)
	figure(fig_id);
	clf;
	plotFMap(C);
	title({ label, sprintf('\\rm\\it(n_L = %.4f, n_D = %.6f)', nL, nD) });
	%plot(cumsum(D .* ((1:size(D,1)).^-1)'));