function [eVec] = massOrthogonal(eVec, M)
  for jj = 1:size(eVec,2) % make sure the eigenvectors are mass-orthogonal
     eVec(:,jj) =  eVec(:,jj) / sqrt((eVec(:,jj)'*M*eVec(:,jj)));
  end
end

