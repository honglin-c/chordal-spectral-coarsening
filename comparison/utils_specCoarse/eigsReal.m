function [eVal, eVec] = eigsReal(L, massMat, numEigs)
[eVec, eVal] = eigs(L + 1e-8.*speye(size(L,1),size(L,1)), massMat, numEigs, 'sm');
[~ ,i] = sort(diag(abs(eVal))); % sort
eVal = sum(eVal,2);
eVal = eVal(i) - 1e-8;
eVec = eVec(:,i);