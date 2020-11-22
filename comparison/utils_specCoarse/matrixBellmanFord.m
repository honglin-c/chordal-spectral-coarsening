function [dist, nearestCenter] = matrixBellmanFord(SoC, centersIdx)
% MATRIXBELLMANFORD performs matrix bellman ford algorithm to compute the
% shortest graph distance for every node to a given set of center points
%
% Reference: 
% Bell, Algebraic Multigrid for Discrete Differential Forms, 2008
%
% [dist, nearestCenter] = matrixBellmanFord(A, centerIdx)
%
% Inputs:
%   SoC |V| x |V| matrix of strength of connections 
%   centersIdx #centers vector of indices of the centers
% Outputs:
%   dist a |V| vector of the graph distance to the closest center point
%   nearestCenter a |V| vector of indices of the nearest center point


% initialize dist and nearestCenter
nV = size(SoC,1);
dist = Inf(nV,1);
nearestCenter = zeros(nV,1);

% set center indices
dist(centersIdx) = 0;
nearestCenter(centersIdx) = centersIdx;

% iterate through all the non-zeros in the system matrix A
[i,j,dij] = find(SoC);

while true
    idx = find( (dist(i)+dij) < dist(j));
    if size(idx,1) ~= 0
        dist(j(idx)) = dist(i(idx)) + dij(idx);
        nearestCenter(j(idx)) = nearestCenter(i(idx));
    else
        break;
    end 
end
end