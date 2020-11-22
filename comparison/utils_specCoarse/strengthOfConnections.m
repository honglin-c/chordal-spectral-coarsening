function SoC = strengthOfConnections(L, M, p)
% STRENGTHOFCONNECTIONS computs the strength of connections (SOC) between
% adjacent vertices given the input operator L and the mass matrix M
% 
% SoC = strengthOfConnections(L,M)
%
% Inputs:
%   L m x m PSD matrix of a differential operator 
% .   (diagonal positive, off-diagonal mostly negative)
%   M m x m diagonal matrix of variable masses (e.g. vertex areas)
% . p a real number controls the unit (e.g. 0.5 for LB operator)
% Outputs:
%   SoC a m x m matrix where off-diagonals are the SoC between variables

if nargin < 3
    p = 1/2;
end

MList = diag(M);
L(1:size(M,1)+1:end) = 0; % remove diagonal
[i,j,val] = find(-L);

% set negative off-diag entries to zeros
negIdx = find(val < 0);
i(negIdx) = [];
j(negIdx) = [];
val(negIdx) = [];

strength = max((MList(i) + MList(j)).^(p) ./ val, 0);
SoC = sparse(i,j,strength, size(L,1), size(L,1), length(val));
SoC = (SoC + SoC') / 2;