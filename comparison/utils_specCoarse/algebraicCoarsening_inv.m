function [Lc, Mc, G, P, K, Cpt, errorHis, time] = algebraicCoarsening_inv(L, M, numNc, varargin)
%%
% default parameter values
lr = 1e-2; % learning rate (step size) for the gradient descent
lrDecayIter = 3; % learning rate decay iterations (optional)
maxIter = 1000; % the maximum iteration (usually won't reach)
stallIter = 5; % stall iteration
landmarks = []; % vertices you want to keep in the coarsening
numEig = round(numNc * 0.5); % number of eigenvectors in use 

% process user input parameters
params_to_variables = containers.Map(...
    {'lr', 'decayIter', 'stallIter', 'maxIter', 'landmarks','numEig'}, ...
    {'lr', 'lrDecayIter', 'stallIter', 'maxIter', 'landmarks','numEig'});
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

%% Combinatorial Coarsening
invM = diag(diag(M).^-1);
SoC = strengthOfConnections(L, M);

% clustering
shuffle = @(v)v(randperm(numel(v)));
seedsIdx = [1:size(L,1)];
seedsIdx(landmarks) = [];
seedsIdx = shuffle(seedsIdx);
seedsIdx = [landmarks, seedsIdx];
seedsIdx = seedsIdx(1:numNc);
clusterAssignment = graphKMediods(SoC, seedsIdx, landmarks);

% construct P (projection) and K (assignment)
[Cpt, ~, idxToCpt] = unique(clusterAssignment);
P = sparse([1:length(Cpt)], Cpt, ones(length(Cpt),1), ...
    length(Cpt), size(L,1), length(Cpt));
K = sparse(idxToCpt, [1:size(L,1)], ones(length(idxToCpt),1), ...
    length(Cpt), size(L,1), length(idxToCpt));

%% Operator Optimization
Mc = K * M * K'; % coarse mass matrix
invMc = diag(diag(Mc).^-1);

S = double(L~=0);
J = double(K * S * K' ~= 0);
H = double(K' * J ~= 0);
A = double(H' * S * H ~= 0);

% eigen vector as test vectors
[eVal, U] = eigsReal(L, M, numEig+1);
eVal_scale = [ones(100,1); 0.01*ones(length(eVal)-1-100,1)];
eVal_inv = diag(1./eVal(2:end).*eVal_scale);
% eVal_inv = diag(1./eVal(2:end));
U = U(:,2:end);

% G = K'; % initial G
G = P'; % initial G use R seems to be better???

% energy function/gradient precomputation
A = P * invM * L * U * eVal_inv;
B = P * U * eVal_inv;
% A = P * invM * L * U;
% B = P * U;
BB = B*B';
AB = A*B';

% projection precomputation (G=G+(U0-G*PU)*invUPPU*PU')
% U0 = U(:,1);
% PU = P*U0;
% invUPPU = inv(U0'*P'*P*U0);
U0 = U(:,1);
ZIdx = find(H);
[rIdx, cIdx] = find(H);
Z = sparse(ZIdx, [1:length(ZIdx)], ones(length(ZIdx),1), size(G,1)*size(G,2), length(ZIdx));
g = G(ZIdx);
Aproj = kron((P*U0)',speye(size(G,1))) * Z;
AAA = Aproj'*inv(Aproj*Aproj');

% stopping criteria (stall)
objValOld = Inf;
stop = 0;
lrDecayCount = 0;
decayRatio = 0.2;

% NADAM
timeStep = 0;
beta1 = 0.9;
beta2 = 0.8;
eps = 1e-8;
mt = sparse(size(G,1), size(G,2));
nt = sparse(size(G,1), size(G,2));

% function for computing objective function
objFunc = @(params) objFuncTemplate(params,A, B, L, Mc, invMc);
errorHis = [];

fprintf('iter 0, cost %.4e\n', objFunc(G));

tic;
for iter = 1:maxIter
    % precomputation
    objVal = objFunc(G);
    errorHis = [errorHis objVal];

    % compute gradient 
    % Note that MATLAB implementation DOES NOT exploit the sparse gradient
    % because looping over the indices is slow in MATLAB. Therefore, the 
    % MATLAB gradient computation is much slower
    LG = L*G;
    BBGLGinvMc = BB*(G'*(LG*invMc));
    grad = LG * (- AB' - AB + BBGLGinvMc + BBGLGinvMc') ;
    grad = grad .* H;

    % NADAM
    timeStep = timeStep + 1;
    mt = beta1*mt + (1-beta1)*grad;
    nt = beta2*nt + (1-beta2)*(grad.^2);
    mt_hat = mt / (1-beta1^timeStep);
    nt_hat = nt / (1-beta2^timeStep);
    Ngrad = 1./(sqrt(nt_hat) + eps) .* (beta1*mt_hat + (1-beta1)*grad/(1-beta1^timeStep));

    % update G
    G = G - lr * Ngrad;
%     G = projFunc(G);
    g1 = G(ZIdx);
    g = g1 - AAA * (Aproj*g1-U0(:));
    G = sparse(rIdx, cIdx, g, size(H,1), size(H,2));
    
%     fprintf('%.1e\n', max(abs(G * P * U0 - U0)))

    % print progress
    if mod(iter,10) == 0
        fprintf('iter %i, cost %.4e\n', iter, objVal);
    end

    % stopping criteria (if the energy stop decreasing for a few iterations)
    if objValOld < objVal
        stop = stop + 1;
        if stop > stallIter
            if lrDecayCount < lrDecayIter % not decrease lr yet
                lr = lr * decayRatio;
                stop = 0;
                lrDecayCount = lrDecayCount + 1;
                fprintf('derease learning rate to: %f\n', lr);
            else
                fprintf('iter %i, cost %.4e\n', iter, objValOld);
                break;
            end
        end
    else
        stop = 0;
        objValOld = objVal;
        Gbest = G;
    end
end
time = toc;
G = Gbest;
Lc = G' * L * G;
Lc = (Lc + Lc')/2; % sometimes a slightly non-symmetric matrix fails 'eigs' due to numerical issues