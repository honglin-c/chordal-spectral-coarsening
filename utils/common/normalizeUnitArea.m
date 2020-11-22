function V = normalizeUnitArea(V,F)
% NORMALIZEUNITAREA normalizes a given mesh to have unit total surface area
%
% V = normalizeUnitArea(V,F)
%
% Inputs:
%   V |V| x 3 matrix of vertex positions
%   F |F| x 3 matrix of indices of triangle corners
% Outputs:
%   V a new matrix of vertex positions whose total area is 1

VA = vertexAreas(V,F);
totalArea = sum(VA);
V = V / sqrt(totalArea); % normalize shape to have unit area