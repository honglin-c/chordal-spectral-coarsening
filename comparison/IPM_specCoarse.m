function [X, out] = IPM_specCoarse(solver, Af, Bf, Cf, D)
    % IPM_specCoarse : Interior point method wrapper (using cvx) for projection onto sparse
    % matrix cones.
    %
    %   [X, out] = IPM((solver,problemtype, A, D): projects sparse matrix 
    %   A onto matrix cones using CVX (interior point solver). Returns X 
    %   the projected matrix and out, a structure containing runtime details.
    %
    % INPUTS
    %   solver      : string indicating which cvx solver to use 
    %                 (e.g. sedumi, sdpt3, mosek). Solver must already be 
    %                 installed.
    %   D           : Sparsity pattern adjacency matrix
    % OUTPUTS
    %   X   : Projection solution
    %   out : Tracked details. Fields include
    %       runtime : CPU time used in solving. Does not include
    %                 preprocessing.
    %       obj     : Final objective value.
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015
    
    N = size(D,1);
    vec1 = ones(N,1);
    vec0 = zeros(N,1);
    cvx_solver(solver)
    
    % solve (||Af - Bf*X*Cf||^2)/2
    % s.t X \in PSD and S_D
    %     X * ones(N,1) = zeros(N,1)
    tic;
    cvx_begin 
      cvx_precision high
      variable X(N,N) symmetric
%       minimize (1/2*sum_square(vec((X-A).*D)))
      minimize (1/2*sum_square(vec(Af - Bf * (X.*D) * Cf)))
      subject to
          X == semidefinite(N)
          X.*D == X
          X * vec1 == vec0
    cvx_end
    endtime = toc;

    out.runtime = endtime;
    out.obj = cvx_optval;
    out.iter = cvx_slvitr;
    out.tol = cvx_slvtol
end
