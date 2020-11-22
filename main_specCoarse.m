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
% objname = './data/models/cell.obj';
objname = './data/models/spot.obj';
% objname = './data/models/monster_110K.ply';
% objname = './data/models/cactus_25K.ply';
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

%% chordal decomposition
Xsp = (Lc~=0);
nV = size(Xsp,1);
X_init = R*L*R'; % initialize the optimization variable to be the same as [Liu et al. 2019]
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
% a = Af(:);    % the lines commented out are another version of the implementation different from Appendix E 
% E = kron(Cf',Bf) * P.x2X;
% G = kron(ones(1,nV), speye(nV)) * P.x2X; % G*x = X*1 = 0
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
t_admm = tic;
[xy,z,state] = admmSolve(argminX,argminZ,admm_A,admm_B,admm_c,state,energyFunc,updatedLHS);
time_admm = toc(t_admm);
x = xy(1:nx);

%% check solution
X = reshape(P.x2X*x, nV,nV);
matrixChecks(X,Mc,'X');

oursVal = sum(sum((Af - Bf * X * Cf).^2)) / 2;
fprintf("our energy: %.4e\n", oursVal);
fprintf("our solve time: %.4e, [Liu et al. 2019] solve time: %.4e\n", time_admm, time_other);  

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