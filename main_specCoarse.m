clear all; close all;
addpath(genpath('./gptoolbox'))
addpath(genpath('./utils'))
addpath(genpath('./comparison'))
rng(0)

% coarsening parameter
numVc = 600;
numEig = round(numVc * 0.5);
% numEig = 100;
numEig_vis = 100; % number of eigs to visualize
normalizationScale = sqrt(numVc);

%% read mesh & gen spectral coarsening data
% objname = './data/models/armadillo.obj';
objname = './data/models/cell.obj';
% objname = './data/models/monster_110K.ply';
% objname = './data/models/cactus_25K.ply';
% objname = './data/models/spot.obj';
[filepath,name,ext] = fileparts(objname);

if strcmp(ext, '.obj')
  [V,F] = readOBJ(objname);
elseif strcmp(ext, '.stl')
  [V,F] = readSTL(objname);
elseif strcmp(ext, '.ply')
  [V,F] = readPLY(objname);
end
V = normalizationScale * normalizeUnitArea(V,F);
L = -cotmatrix(V,F);
M = massmatrix(V,F);
tic;
[eVal,eVec] = eigsReal(L, M, numEig+1);
toc;

% run [Liu et al. 2019] to generate sparsity pattern
fprintf("#high-res mesh vertices = %d, #coarse mesh vertices = %d\n", size(V,1), numVc);
fprintf('===============[Liu et al. 2019] starts===============\n');
[Lc, Mc, ~, R, K, Cpt, errorHis, time_other] = algebraicCoarsening(L, M, numVc, 'numEig', numEig,'maxIter', 1000);
fprintf('===============[Liu et al. 2019] finishes===============\n');

% set the sparsity pattern
Xsp = (Lc~=0);
X_init = R*L*R'; % initialize the optimization variable to be the same as [Liu et al. 2019]

%% ADMM solve with chordal decomposition
[X, state] = chordal_coarsening(L, M, R, Xsp, Mc, X_init, 'numEig', numEig);

fprintf("our energy: %.4e\n", state.energy_ours);
fprintf("our solve time: %.4e, [Liu et al. 2019] solve time: %.4e\n", state.time_admm, time_other);  

%% visualize the results
% plot results
[eVal,eVec] = eigsReal(L, M, numEig_vis);
[eVal_ours, eVec_ours] = eigsReal(X, Mc, numEig_vis);
[eVal_spec, eVec_spec] = eigsReal(Lc, Mc, numEig_vis);
for jj = 1:numEig_vis % make sure eigenvectors are orthogonal
   eVec_ours(:,jj) =  eVec_ours(:,jj) / sqrt((eVec_ours(:,jj)'*Mc*eVec_ours(:,jj)));
   eVec(:,jj) =  eVec(:,jj) / sqrt((eVec(:,jj)'*M*eVec(:,jj)));
   eVec_spec(:,jj) =  eVec_spec(:,jj) / sqrt((eVec_spec(:,jj)'*Mc*eVec_spec(:,jj)));
end
fMap_ours = eVec_ours' * Mc * R * eVec;
fMap_spec = eVec_spec' * Mc * R * eVec;
figure(1)
subplot(1,2,1)
plotFMap(fMap_ours)
title('fmap (ours)')
subplot(1,2,2)
plotFMap(fMap_spec)
title('fmap [Liu et al. 2019]')

figure(2)
plot(log10(state.objHis))
hold on
plot(log10(errorHis))
xlabel('iterations')
ylabel('log10 energy value')
legend('ours','[Liu et al. 2019]')

figure(3)
plot(eVal)
hold on
plot(eVal_ours)
plot(eVal_spec)
legend('original', 'ours','[Liu et al. 2019]')
ylabel('eigenvalues')

% visualize one eigenfunctions (this may have sign flip)
figure(4)
subplot(1,3,1)
plotMesh(V,F,eVec(:,30))
title('original');
subplot(1,3,2)
scatter3(V(Cpt,1),V(Cpt,2),V(Cpt,3), 30, eVec_ours(:,30), 'filled')
title('ours');
axis equal off
subplot(1,3,3)
scatter3(V(Cpt,1),V(Cpt,2),V(Cpt,3), 30, eVec_spec(:,30), 'filled')
title('[Liu et al. 2019]');
axis equal off
sgtitle('visualize one eigen function (may have sign flip)')

figure(5)
subplot(1,2,1)
spy(Lc);
title('[Liu et al. 2019]');
subplot(1,2,2)
spy(X)
title('our pattern');

eVecc = R * eVec;
E_L = sum(sum(eVec' * L * eVec));
E_X = sum(sum(eVecc' * X * eVecc));
E_Lc = sum(sum(eVecc' * Lc * eVecc));
fprintf('E_L = %.4e, E_X = %.4e, E_Lc = %.4e\n', E_L, E_X, E_Lc);