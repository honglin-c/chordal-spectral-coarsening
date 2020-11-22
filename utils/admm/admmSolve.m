function [X,Z,state] = admmSolve(argmin_X,argmin_Z,A,B,c,state,energyFunc,updatedLHS)
    % a modified admm.m from gptoolbox
    % https://github.com/alecjacobson/gptoolbox/blob/master/matrix/admm.m
    
    if(~exist('energyFunc','var'))
        energyFunc = [];
    else
        objHis = [];
    end
    
    tol_abs = 5e-6;
    tol_rel = 1e-4;
    check_interval = 5;
    print_interval = 10;
    bmu = 5;
    btao_inc = 2;
    btao_dec = 2;
    max_iter = 300;

    % Initial conditions
    if ~isfield(state, 'U')
        if ~isfield(state, 'X')
            state.X = rand(size(A,2),size(c,2));
        end
        if ~isfield(state, 'Z')
            state.Z = rand(size(B,2),size(c,2));
        end
        state.U = zeros(size(c,1),size(c,2));
        state.rho_prev = nan;
        state.rho = 2e-2;
    end
    
    % cache the decomposition
    state.dLHS = updatedLHS(state.rho);
    
    % print the energy of 0-th iteration before the optimization
    nz = length(state.Z);
    objHis = [objHis; full(energyFunc(state.X(1:end-nz)))];
    fprintf('iter: 0, obj: %.4e \n', objHis(end));
    
    state.argmin_x_time = 0;
    state.argmin_z_time = 0;

    cnorm = norm(c,'fro');
    for iter = 1:max_iter    
        tx = tic;
        [state.X] = argmin_X(state.Z,state.U,state.rho,state.dLHS);
        state.argmin_x_time = state.argmin_x_time + toc(tx);
        
        state.Z_prev = state.Z;
        tz = tic;
        [state.Z] = argmin_Z(state.X,state.U,state.rho);
        state.argmin_z_time = state.argmin_z_time + toc(tz);
        
        state.U_prev = state.U;
        state.U = state.U+A*state.X+B*state.Z-c;
        state.rho_prev = state.rho;

        dual_residual = state.rho*norm(A'*B*(state.Z_prev - state.Z),'fro');
        residual = norm(A*state.X+B*state.Z-c,'fro');
        if mod(iter,check_interval) == 0
            if residual > bmu*dual_residual
                state.rho = btao_inc*state.rho;
                % From python code: https://github.com/tneumann/cmm/blob/master/cmmlib/cmm.py
                state.U = state.U/btao_inc;
                % cache the decomposition
                state.dLHS = updatedLHS(state.rho);
            elseif dual_residual > bmu*residual
                state.rho = state.rho/btao_dec;
                state.U = state.U*btao_dec;
                % cache the decomposition
                state.dLHS = updatedLHS(state.rho);
            end
        end
%         % From Python code (this seems weird)
%         k = size(c,2);
%         eps_pri = sqrt(k*2)*tol_abs + tol_rel*max([norm(A*state.X,'fro'),norm(B*state.Z,'fro'),cnorm]);
%         eps_dual = sqrt(k)*tol_abs +  tol_rel*state.rho*norm(A'*state.U,'fro');
        
        p = size(A,1);
        n = size(A,2);
        eps_pri  = sqrt(p)*tol_abs + tol_rel*max([norm(A*state.X,'fro'), norm(B*state.Z,'fro'), cnorm]);
        eps_dual = sqrt(n)*tol_abs + tol_rel*norm(A'*state.U,'fro');
        
        if residual < eps_pri && dual_residual < eps_dual
            break;
        end
        
        if isempty(energyFunc)
            fprintf('iter:%d, r_pri: %.4e, r_dual: %.4e, rho: %.4e \n', iter, residual, dual_residual, state.rho)
        else
            nz = length(state.Z);
            objHis = [objHis; full(energyFunc(state.X(1:end-nz)))];
            if mod(iter,print_interval) == 0
              fprintf('iter:%d, obj: %.4e, r_pri: %.4e, r_dual: %.4e, rho: %.4e \n', iter, objHis(end), residual, dual_residual, state.rho)
            end
        end
    end
    X = state.X;
    Z = state.Z;
    if ~isempty(energyFunc)
        state.objHis = objHis;
    end
end