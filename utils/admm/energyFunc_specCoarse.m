function energy = energyFunc_specCoarse(x,Af,Bf,Cf,x2X,nV)
  X = reshape(x2X*x,nV,nV);
  energy = sum(sum((Af - Bf * X * Cf).^2)) / 2;
end