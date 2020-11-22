function VA = vertexAreas(V, F)
% VERTEXAREAS computs per vertex area of a triangle mesh
% 
% VA = vertexAreas(V,F)
%
% Inputs:
%   V |V| x 3 matrix of vertex positions
%   F |F| x 3 matrix of indices of triangle corners
% Outputs:
%   VA a |V| vector of vertex areas (summation of 1/3 adjacent face areas)

FN = cross(V(F(:,1),:)-V(F(:,2),:), V(F(:,1),:) - V(F(:,3),:));
FA = sqrt(sum(FN.^2,2)) ./ 2; % face area

rIdx = [F(:,1);F(:,2);F(:,3)];
cIdx = ones(size(rIdx));
val = [FA(:,1);FA(:,1);FA(:,1)] ./ 3;

VA = full(sparse(rIdx,cIdx,val,size(V,1),1));