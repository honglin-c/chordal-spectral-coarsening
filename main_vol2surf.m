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

%% chordal decomposition
% get sparsity
Xsp = (Lc~=0);
Xsp = Xsp^numRing;

%% chordal decomposition
nV = size(Xsp,1);
X_init = Ls; 
tfill = 200; % parameter for chordal embedding
tsize = 200; % parameter for chordal embedding
[X_chordal, cliques, vecCliques, P] = preprocessChordal(Xsp, tfill, tsize);

%% vectorize spectral coarsening energy 
% (matrix version) energy = 0.5 * || Af - Bf * X * Cf ||^2_F
% (sparse vector version) energy = 0.5 * || a - E * x ||^2_F
% where x is the lower triangular non-zeros of X (assume X is symmetric)
[eVal,eVec] = eigsReal(L, M, numEig+1);
eVec = massOrthogonal(eVec, M); % make sure eigenvectors are orthogonal
eVec = eVec(:,2:end);
eVal_scale = [ones(100,1); 0.1*ones(length(eVal)-1-100,1)];
eVal_inv = diag(1./eVal(2:end).*eVal_scale);
% eVal_inv = diag(1./eVal(2:end));
invM = diag(diag(M).^(-1));
invMc = diag(diag(Mc).^(-1));
sqrtMc = diag(sqrt(diag(Mc)));
Af = sqrtMc * R * invM * L * eVec * eVal_inv;
Bf = sqrtMc * invMc;
Cf = R * eVec * eVal_inv;
% Af = sqrtMc * R * invM * L * eVec;
% Bf = sqrtMc * invMc;
% Cf = R * eVec;
% a = Af(:);
% E = kron(Cf',Bf) * P.x2X;
EE = P.x2X' * kron(Cf*Cf',Bf'*Bf) * P.x2X;
Ea_tmp = Bf'* Af * Cf';
Ea = P.x2X' * Ea_tmp(:);
G = kron(ones(1,nV), speye(nV)) * P.x2X; % G*x = X*1 = 0

%% prepare ADMM
nx = size(P.X2x,1);
nxc = size(P.x2xc,1);
nz = size(P.z2xc,2);
nG = size(G,1);

% get the argmin data
argmin_data.nV = nV;
argmin_data.nx = nx;
argmin_data.nxc = nxc;
argmin_data.nz = nz;
argmin_data.nG = nG;

argmin_data.G = G;
argmin_data.Ea = Ea;
argmin_data.Mc = Mc;
argmin_data.P = P;
argmin_data.cliques = cliques;
argmin_data.vecCliques = vecCliques;

LHS_fixed = [EE, P.x2xc', G'; 
             P.x2xc, sparse(nxc, nxc), sparse(nxc, nG);
             G, sparse(nG,nxc), sparse(nG, nG)];
[row,col,value] = find(P.z2xc * P.z2xc');
LHS_PxcPxc = sparse(row+nx, col+nx, value, ...
                    size(LHS_fixed,1), size(LHS_fixed,2));
                  
fprintf('LHS equation number = %d\n', size(LHS_fixed,1)); 

% get the argmin functions
argminX = @(z,u,rho,dLHS) argminX_template(z,u,rho,dLHS,argmin_data);
argminZ = @(xy,u,rho) argminZ_template(xy,u,rho,argmin_data);
updatedLHS = @(rho) update_argminX_dLHS(LHS_fixed, LHS_PxcPxc, rho);
energyFunc = @(x) energyFunc_specCoarse(x,Af,Bf,Cf,P.x2X,nV);

% min f(x) + g(z)
% st Ax + Bz + c = 0
admm_A = [sparse(nz,nx), speye(nz,nz)];
admm_B = -speye(nz,nz);
admm_c = sparse(nz,1);

% initialize ADMM
state.Z = P.X2z_avg * X_init(:);
state.X = [P.X2x * X_init(:); state.Z];

%% admm solve
tic;
[xy,z,state] = admmSolve(argminX,argminZ,admm_A,admm_B,admm_c,state,energyFunc,updatedLHS);
toc;
x = xy(1:nx);

%% check solution
X = reshape(P.x2X*x, nV,nV);
matrixChecks(X,Mc,'X');

cotanVal = sum(sum((Af - Bf * Lc * Cf).^2)) / 2;
oursVal = sum(sum((Af - Bf * X * Cf).^2)) / 2;
fprintf("our energy: %.4e, cotan energy: %.4e\n", oursVal, cotanVal);

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