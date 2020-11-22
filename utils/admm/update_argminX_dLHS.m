function dLHS = update_argminX_dLHS(LHS_fixed, LHS_rho, rho)
  LHS = LHS_fixed - 1/rho * LHS_rho;
  dLHS = decomposition(LHS, 'ldl');
  flag = isIllConditioned(dLHS);
  if flag
    warning('The decomposition is ill-conditioned');
  end
end