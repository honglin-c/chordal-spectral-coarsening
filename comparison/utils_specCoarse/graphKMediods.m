function [clusterAssignIdx] = graphKMediods(SoC, seedsIdx, landmarks, maxIter)
% GRAPHKMEDIODS performs k-mediods clustering given a strength of
% connection matrix (see strengthOfConnections.m)
% 
% clusterAssignIdx = graphKMediods(SoC, seedsIdx)

% Inputs:
%   SoC |V| x |V| matrix of strength of connections 
%   seedsIdx #seeds vector of indices of initial seeds for the clustering
%   maxIter maximum iterations to perform coarsening
% Outputs:
%   clusterAssignIdx a |V| vector of cluster assignment result

tStart = tic;

if (nargin < 3)
    landmarks = [];
end
if (nargin < 4)
    maxIter = 20;
end

clusterAssignIdx = zeros(size(SoC,1),1);
[row, col, dij] = find(SoC);
ijList = unique(sort([row, col],2),'rows');
for iter = 1:maxIter
    [~, nearestCenter] = matrixBellmanFord_fast(SoC, seedsIdx, row, col, dij);
    isBorder = find(nearestCenter(ijList(:,1)) ~= nearestCenter(ijList(:,2)));
    borderIdx = ijList(isBorder,:);
    borderIdx = unique(borderIdx(:));

    [distFromBorder, ~] = matrixBellmanFord_fast(SoC, borderIdx, row, col, dij);
    for ii = (length(landmarks)+1):length(seedsIdx)
        cIdx = seedsIdx(ii);
        aggNodes = find(nearestCenter == cIdx);
        [~,newCIdx] = max(distFromBorder(aggNodes));
        seedsIdx(ii) = aggNodes(newCIdx);
    end
    
    if clusterAssignIdx == nearestCenter
        break;
    else
        clusterAssignIdx = nearestCenter;
    end
end

fprintf('graphKMediods iteration %d\n', iter)
tEnd = toc(tStart);
fprintf('graphKMediods %.4f sec\n', tEnd)
end

function [dist, nearestCenter] = matrixBellmanFord_fast(A, centerIdx, i, j, dij)
% MATRIXBELLMANFORD_FAST performs matrix bellman ford algorithm to compute
% graph distance for every node to a given set of seed points
%
% Note:
% this is a simplified version designed specifically for the graphKMediods.
% For a complete version, please see matrixBellmanFord.m
%
% Reference: 
% Bell, Algebraic Multigrid for Discrete Differential Forms, 2008

% initialize dist and nearestCenter
nV = size(A,1);
dist = Inf(nV,1);
nearestCenter = zeros(nV,1);

% set center indices
dist(centerIdx) = 0;
nearestCenter(centerIdx) = centerIdx;

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