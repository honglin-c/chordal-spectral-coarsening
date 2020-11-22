clear all; close all;
addpath(genpath('./gptoolbox'))
addpath(genpath('./utils'))
addpath(genpath('./comparison'))

% coarsening parameter
numVc = 500;
numEig = 100;
numEig_vis = 100; % number of eigs to visualize
normalizationScale = sqrt(numVc);
numRings = 1; % number of rings (change this to set 1-, 2- and 3-ring)

%% read mesh
% objname = './data/models/armadillo.obj';
objname = './data/models/spot.obj';
[filepath,name,ext] = fileparts(objname);

if strcmp(ext, '.obj')
  [V,F] = readOBJ(objname);
elseif strcmp(ext, '.stl')
  [V,F] = readSTL(objname);
elseif strcmp(ext, '.ply')
  [V,F] = readPLY(objname);
end
V = normalizationScale * normalizeUnitArea(V,F);

%% generate coarse mesh using [Garland & Heckbert 1997]
[Vc,Fc,R] = QSlim_withR(V,F,numVc);
fprintf('============Finish generating coarse mesh connection using [Garland & Heckbert 1997]==========\n');
fprintf('==================================Our algorithm starts========================================\n');

%% compute cotangent laplacian and mass matrix
L = -cotmatrix(V,F);
M = massmatrix(V,F);
Lc = -cotmatrix(Vc,Fc);
Mc = massmatrix(Vc,Fc);

% set the sparsity pattern
Xsp = (Lc~=0);
Xsp = (Xsp^numRings ~= 0); % sparsity pattern
X_init = Lc; 

%% ADMM solve with chordal decomposition
[X, state] = chordal_coarsening(L, M, R, Xsp, Mc, X_init, 'numEig', numEig);

fprintf("our solve time: %.4e\n", state.time_admm);  

%% visualize the results
% plot results
[eVal,eVec] = eigsReal(L, M, numEig);
[eVal_ours, eVec_ours] = eigsReal(X, Mc, numEig);
[eVal_cotan, eVec_cotan] = eigsReal(Lc, Mc, numEig);
% make sure the eigenvectors are mass-orthogonal
eVec_ours = massOrthogonal(eVec_ours, Mc);
eVec = massOrthogonal(eVec, M);
eVec_cotan = massOrthogonal(eVec_cotan, Mc);
fMap_ours = eVec_ours' * Mc * R * eVec;
fMap_cotan = eVec_cotan' * Mc * R * eVec;
figure(1)
subplot(1,2,1)
plotFMap(fMap_ours(1:numEig_vis, 1:numEig_vis))
title('fmap (ours)')
subplot(1,2,2)
plotFMap(fMap_cotan(1:numEig_vis, 1:numEig_vis))
title('fmap cotan Laplace')

figure(2)
plot(log10(state.objHis))
xlabel('iterations')
ylabel('log10 energy value')

figure(3)
plot(eVal)
hold on
plot(eVal_ours)
plot(eVal_cotan)
legend('original', 'ours','cotan Laplace')
ylabel('eigenvalues')

% visualize one eigenfunctions (this may have sign flip)
figure(4)
subplot(1,3,1)
plotMesh(V,F,eVec(:,40))
title('original');
subplot(1,3,2)
plotMesh(Vc,Fc,eVec_ours(:,40))
title('ours');
subplot(1,3,3)
plotMesh(Vc,Fc,eVec_cotan(:,40))
title('[Garland & Heckbert 1997]');
axis equal off
sgtitle('visualize one eigen function (may have sign flip)')

figure(5)
subplot(1,2,1)
spy(Lc);
title('[Garland & Heckbert 1997]');
subplot(1,2,2)
spy(X)
title('our pattern');

% visualize functional map
fmap = fmap_operator(numEig_vis, L, M);

other = struct();
[other.C, other.nC, other.nL, other.nD] = fmap(R, Lc, Mc);
fig_fmap(6, '[Garland & Heckbert 1997]', other.C, other.nL, other.nD);

our = struct();
[our.C, our.nC, our.nL, our.nD] = fmap(R, X, Mc);
fig_fmap(7, '[Ours]', our.C, our.nL, our.nD);
