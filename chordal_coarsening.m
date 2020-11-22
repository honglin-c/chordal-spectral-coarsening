function [X, state] = chordal_coarsening(L, M, R, Xsp, Mc, X_init, varargin)
% CHORDAL_COARSENING 
% coarsen a sparse positive semidefinite discrete operator from #V by #V 
% to #Vc by #Vc while preserving its spectral properties, i.e., the
% eigenvectors and eigenvalues at the lowes
%
% Inputs:
%   L       #V by #V    sparse, positive semi-definite high-res discrete operator (e.g., laplacian matrix) 
%   M       #V by #V    mass matrix of the high-res operator
%   R       #Vc by #V   restriction operator
%   Xsp     #Vc by #Vc  sparsity pattern of the coarsened operator
%   Mc      #Vc by #Vc  mass matrix of the coarsened operator
%   X_init  #Vc by #Vc  initial guess of the coarsened operator
%   Optional:
%     'numEig'  number of the eigenvectors in use
%     'tfill'   parameter controlling the chordal decomposition
%     'tsize'   (another) parameter controlling the chordal decomposition
% Outputs:
%   X       #Vc by #Vc  the output coarsened operator after optimization    
%   state   optimization log (optional)
%
% See: % "Chordal decomposition for spectral coarsening" [Chen et al. 2020]
%      https://github.com/honglin-c/chordal-spectral-coarsening
%

  nV = size(Xsp, 1);
  numVc = size(Xsp, 1);
  numEig = round(numVc * 0.5);

  % process user input parameters
  params_to_variables = containers.Map(...
      {'numEig', 'tfill', 'tsize'}, ...
      {'numEig', 'tfill', 'tsize'});
  v = 1;
  while v <= numel(varargin)
      param_name = varargin{v};
      if isKey(params_to_variables,param_name)
          assert(v+1<=numel(varargin));
          v = v+1;
          % Trick: use feval on anonymous function to use assignin to this workspace 
          feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
      else
          error('Unsupported parameter: %s',varargin{v});
      end
      v=v+1;
  end

  %% perform chordal decomposition
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

  if numEig <= numVc
    a = Af(:); 
    E = kron(Cf',Bf) * P.x2X;
    G = kron(ones(1,numVc), speye(numVc)) * P.x2X; % G*x = X*1 = 0
    EE = E'*E;
    Ea = E'*a;
  else  % if too many eigenvectors are in use, choose another version of the vectorization (see Appendix E)
    G = kron(ones(1,nV), speye(nV)) * P.x2X; % G*x = X*1 = 0
    EE = P.x2X' * kron(Cf*Cf',Bf'*Bf) * P.x2X;
    Ea_tmp = Bf'* Af * Cf';
    Ea = P.x2X' * Ea_tmp(:);
  end

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
  t_admm = tic;
  [xy,z,state] = admmSolve(argminX,argminZ,admm_A,admm_B,admm_c,state,energyFunc,updatedLHS);
  state.time_admm = toc(t_admm);
  x = xy(1:nx);

  %% check solution
  X = reshape(P.x2X*x, nV,nV);
  matrixChecks(X,Mc,'X');

  state.energy_init = sum(sum((Af - Bf * X_init * Cf).^2)) / 2;
  state.energy_ours = sum(sum((Af - Bf * X * Cf).^2)) / 2;

