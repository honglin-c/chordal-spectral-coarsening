clear all; close all;
addpath(genpath('./gptoolbox'))
addpath(genpath('./utils'))
addpath(genpath('./comparison'))

%% read mesh & gen spectral coarsening data
meshname = 'craneHook';
% meshname = 'max-planck';

% read meshes tet mesh (VT,T), surface mesh with links (VF, F)
[VT,T,~] = readMESH(['./data/vol2surf/' meshname '.mesh']);
[VF,F,~] = readMESH(['./data/vol2surf/' meshname '_sf.mesh']);

% coarsening parameter
numVc = size(VF,1);
numEig = size(VT,1)-1;
numRing = 1;
numEig_vis = 100; % number of eigs to visualize
normalizationScale = numVc^(1/3);

% % scale the mesh 
VFpre = VF;
VF = normalizationScale * normalizeUnitVolume(VF,F);
VT = VT * (mean(VF) / mean(VFpre));

% compute laplacian and mass matrix for original tetrahedral mesh
nVT = size(VT,1);
volT = volume(VT,T);
volT_ratio = max(volT) / min(volT)
M = massmatrix(VT,T);
L = -cotmatrix(VT,T);

% compute laplacian for surface mesh with links
nVF = size(VF,1);
volF = volume(VF,F);
volF_ratio = max(volF) / min(volF)
Lc = -cotmatrix(VF,F);

% compute laplacian and mass matrix for original surface mesh
[FT,~,~] = boundary_faces(T); 
[~,IM,VF2VT] = remove_unreferenced(VT,FT);
Fs = IM(FT);
Ls = -cotmatrix(VF,Fs);
Ms = massmatrix(VF,Fs);

Mc = Ms;
Mc = sum(diag(M)) / sum(diag(Mc)) * Mc;

% restriction operator
R = sparse([1:nVF], VF2VT, 1, nVF, nVT);

% get sparsity
Xsp = (Lc~=0);
Xsp = Xsp^numRing;

% set the initial guess of coarsened operator
X_init = Ls; 

%% ADMM solve with chordal decomposition
[X, state] = chordal_coarsening(L, M, R, Xsp, Mc, X_init, 'numEig', numEig, 'weight_scale', true);

fprintf("our solve time: %.4e\n", state.time_admm);  

% visualize the functional map
[eVal,eVec] = eigsReal(L, M, numEig_vis);
[eVal_ours, eVec_ours] = eigsReal(X, Mc, numEig_vis);
[eVal_surf, eVec_surf] = eigsReal(Ls, Mc, numEig_vis);
% make sure eigenvectors are orthogonal
eVec_ours = massOrthogonal(eVec_ours, Mc);
eVec = massOrthogonal(eVec, M); 
eVec_surf = massOrthogonal(eVec_surf, Mc); 
fMap_ours = eVec_ours' * Mc * R * eVec;
fMap_surf = eVec_surf' * Mc * R * eVec;
figure(1)
subplot(1,2,1)
plotFMap(fMap_ours(1:numEig_vis, 1:numEig_vis))
title('fmap (ours)')
subplot(1,2,2)
plotFMap(fMap_surf(1:numEig_vis, 1:numEig_vis))
title('fmap surface cotan Laplace')

% visualize the energy curve
figure(2)
plot(log10(state.objHis))
xlabel('iterations')
ylabel('log10 energy value')

% visualize the eigenvalues
figure(3)
plot(eVal)
hold on
plot(eVal_ours)
plot(eVal_surf)
legend('volumetric', 'ours','surface')
ylabel('eigenvalues')

figure(4)
subplot(1,2,1)
spy(Ls);
title('surface pattern');
subplot(1,2,2)
spy(X)
title('our pattern');