function V = normalizeUnitVolume(V,T)
% NORMALIZEUNITAREA normalizes a given mesh to have unit total volume
%
% V = normalizeUnitVolume(V,T)
%
% Inputs:
%   V |V| x 3 matrix of vertex positions
%   T |T| x 4 matrix of indices of tetrahedron corners
% Outputs:
%   V a new matrix of vertex positions whose total volume is 1

v = volume(V,T);
totalVolume = sum(v);
V = V / totalVolume^(1/3); % normalize shape to have unit volume