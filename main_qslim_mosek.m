clear all; close all;
addpath(genpath('./gptoolbox'))
addpath(genpath('./utils'))
addpath(genpath('./comparison'))

%% The test for mosek requires CVX and MOSEK to run

% coarsening parameter
numVc = 600;  % increase this too much may cause MOSEK to crash your MATLAB because of memory issue
numEig = round(0.5*numVc);
numEig_vis = 100; % number of eigs to visualize
normalizationScale = sqrt(numVc);
numRings = 1;

%% read mesh
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

%% compute cotangent laplacian and mass matrix
L = -cotmatrix(V,F);
M = massmatrix(V,F);
Lc = -cotmatrix(Vc,Fc);
Mc = massmatrix(Vc,Fc);

%% chordal decomposition
Xsp = (Lc~=0);
Xsp = (Xsp^numRings ~= 0); % sparsity pattern
nV = size(Xsp,1);
X_init = Lc; 
tfill = 200; % parameter for chordal embedding
tsize = 200; % parameter for chordal embedding
[X_chordal, cliques, vecCliques, P] = preprocessChordal(Xsp, tfill, tsize);

%% vectorize spectral coarsening energy 
% (matrix version) energy = 0.5 * || Af - Bf * X * Cf ||^2_F
% (sparse vector version) energy = 0.5 * || a - E * x ||^2_F
% where x is the lower triangular non-zeros of X (assume X is symmetric)
[eVal,eVec] = eigsReal(L, M, numEig+1);
eVec = massOrthogonal(eVec, M);
eVec = eVec(:,2:end);
eVal_inv = diag(1./eVal(2:end));
invM = diag(diag(M).^(-1));
invMc = diag(diag(Mc).^(-1));
sqrtMc = diag(sqrt(diag(Mc)));
Af = sqrtMc * R * invM * L * eVec * eVal_inv;
Bf = sqrtMc * invMc;
Cf = R * eVec * eVal_inv;

fprintf('=======================MOSEK starts====================\n');
%% mosek solve
t_mosek = tic;
[Lm, out] = IPM_specCoarse('mosek', Af, Bf, Cf, Xsp);
time_mosek = toc(t_mosek);

fprintf('=======================MOSEK finishes====================\n');
fprintf('====================Our algorithm starts=================\n');

%% admm solve
t_admm = tic;
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
% argmin_data.Ea = E'*a;
argmin_data.Mc = Mc;
argmin_data.P = P;
argmin_data.cliques = cliques;
argmin_data.vecCliques = vecCliques;

LHS_fixed = [EE, P.x2xc', G'; 
             P.x2xc, sparse(nxc, nxc), sparse(nxc, nG);
             G, sparse(nG,nxc), sparse(nG, nG)];
% LHS_fixed = [E'*E, P.x2xc', G'; 
%              P.x2xc, sparse(nxc, nxc), sparse(nxc, nG);
%              G, sparse(nG,nxc), sparse(nG, nG)];
[row,col,value] = find(P.z2xc * P.z2xc');
LHS_PxcPxc = sparse(row+nx, col+nx, value, ...
                    size(LHS_fixed,1), size(LHS_fixed,2));

% get the argmin functions
argminX = @(z,u,rho,dLHS) argminX_template(z,u,rho,dLHS,argmin_data);
argminZ = @(xy,u,rho) argminZ_template(xy,u,rho,argmin_data);
updatedLHS = @(rho) update_argminX_dLHS(LHS_fixed, LHS_PxcPxc, rho);
% energyFunc = @(x) sum((a - E*x).^2) / 2;
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
[xy,z,state] = admmSolve(argminX,argminZ,admm_A,admm_B,admm_c,state,energyFunc,updatedLHS);
time_admm = toc(t_admm);
x = xy(1:nx);

%% check solution
X = reshape(P.x2X*x, nV,nV);
matrixChecks(Lm,Mc,'X');

Lc = Lm;

mosekVal = sum(sum((Af - Bf * Lc * Cf).^2)) / 2;
oursVal = sum(sum((Af - Bf * X * Cf).^2)) / 2;
fprintf("our energy: %.4e, mosek energy: %.4e\n", oursVal, mosekVal);
fprintf("our time: %.4e, mosek time: %.4e\n", time_admm, time_mosek);

% plot results
[eVal,eVec] = eigsReal(L, M, numEig);
[eVal_ours, eVec_ours] = eigsReal(X, Mc, numEig);
[eVal_mosek, eVec_mosek] = eigsReal(Lc, Mc, numEig);
% make sure the eigenvectors are mass-orthogonal
eVec_ours = massOrthogonal(eVec_ours, Mc);
eVec = massOrthogonal(eVec, M);
eVec_mosek = massOrthogonal(eVec_mosek, Mc);
fMap_ours = eVec_ours' * Mc * R * eVec;
fMap_mosek = eVec_mosek' * Mc * R * eVec;
figure(1)
subplot(1,2,1)
plotFMap(fMap_ours(1:numEig_vis, 1:numEig_vis))
title('fmap (ours)')
subplot(1,2,2)
plotFMap(fMap_mosek(1:numEig_vis, 1:numEig_vis))
title('fmap MOSEK')

figure(2)
plot(eVal)
hold on
plot(eVal_ours)
plot(eVal_mosek)
legend('original', 'ours','MOSEK')
ylabel('eigenvalues')

% visualize one eigenfunctions (this may have sign flip)
figure(3)
subplot(1,3,1)
plotMesh(V,F,eVec(:,40))
title('original');
subplot(1,3,2)
plotMesh(Vc,Fc,eVec_ours(:,40))
title('ours');
subplot(1,3,3)
plotMesh(Vc,Fc,eVec_mosek(:,40))
title('MOSEK');
axis equal off

figure(4)
subplot(1,2,1)
spy(Lc);
title('MOSEK pattern');
subplot(1,2,2)
spy(X)
title('our pattern');

% visualize functional map
fmap = fmap_operator(numEig_vis, L, M);

other = struct();
[other.C, other.nC, other.nL, other.nD] = fmap(R, Lc, Mc);
fig_fmap(5, 'MOSEK', other.C, other.nL, other.nD);

our = struct();
[our.C, our.nC, our.nL, our.nD] = fmap(R, X, Mc);
fig_fmap(6, '[Ours]', our.C, our.nL, our.nD);

