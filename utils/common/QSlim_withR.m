function [RV, RF,R] = QSlim_withR(V,F,tarV, varargin)
%     QSlim implement the quadric error edge collapse algorithm [Garland, Heckbert 1997]
%     
%     Inputs:
%     V: |V|-3 vertex list
%     F: |F|-3 face list
%     tarV: target number of simplified vertices
%     
%     Inputs (optional):
%     bdWScale: weights on preserving boundary edges
%     method: "QEM" or "midPoint", using QEM or a simple mid point collapse
%     triQ_threshold: 0~1 triangle quality
%     dotn_threshold: 0~1 normal changes before/after collapse
%     
%     Outputs:
%     RV: |RV|-3 simplified vertices
%     RF: |RF|-3 simplified faces
%
%     Examples:
%     [V,F] = readOBJ('path/to/mesh.obj')
%     [RV, RF] = QSlim(V,F,200);
 
    % parse parameters 
    triQ_threshold = 0.1;
    dotn_threshold = 0.5;
    bdWScale = 3;
    method = 'QEM'; % QEM, midPoint

    params = { ...
        'triangleQuality','NormalDeviation', 'decimationMethod', ...
        'boundaryWeights'};
    vars = { ...
        'triQ_threshold','dotn_threshold', 'method', ...
        'bdWScale'};

    % Map of parameter names to variable names
    params_to_variables = containers.Map(params,vars);
    v = 1;
    while v <= numel(varargin)
        param_name = varargin{v};
        if isKey(params_to_variables,param_name)
            assert(v+1<=numel(varargin));
            v = v+1;
            feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
        else
            error('Unsupported parameter: %s',varargin{v});
        end
            v=v+1;
    end
    
    % check whether input mesh is manifold
    VCheck = is_vertex_nonmanifold(F);
    if any(VCheck)
        error('vertex non-manifold\n');
    end
    [~,ECheck] = nonmanifold_edges(F);
    if any(ECheck)
        error('edge non-manifold\n');
    end
    
    % initialization
    R = speye(size(V,1),size(V,1)); % restriction operator
    [E,ECost,EOptPos,ETag,KV] = edgeQuadric(V,F,bdWScale,method);

    % start edge collapse
    totalCollapse = size(V,1) - tarV;
    numCollapse = 0;
    stall = 0;
    maxStall = 500;
    
    while true
        % print progress
        if mod(numCollapse, 100) == 0
            fprintf('collapse progress [Garland & Heckbert 1997] %d / %d\n', [numCollapse totalCollapse])
        end
        stall = stall + 1;
        
        % get edge with the min cost
        [~,e] = min(ECost);
        isBd = ETag(e); % is boundary or not
        
        % CHECK: if the edge is degenerated
        if E(e,1) == E(e,2) 
            E(e,:) = [];
            ECost(e) = [];
            EOptPos(e,:) = [];
            ETag(e,:) = [];
            stall = 0;
            continue;
        end
        
        % save previous mesh
        Vpre = V;
        Fpre = F;
        
        % get adjacent faces of this edge
        vi = E(e,1);
        vj = E(e,2);
        [adjFi,~] = find(F == vi);
        [adjFj,~] = find(F == vj);
        adjF = unique([adjFi; adjFj]);
        
        % CHECK: check link condition
        Nvi = unique(F(adjFi,:));
        Nvj = unique(F(adjFj,:));
        if (~isBd) && (length(intersect(Nvi, Nvj))~=4) % 4 includes vi, vj
            fprintf('(interior) link condition fail \n');
            ECost(e) = inf;
            if stall > maxStall break; end
            continue;
        end
        
        % compute face normals before edge collapse
        FN_prev = cross(V(F(adjF,1),:)-V(F(adjF,2),:), V(F(adjF,1),:) - V(F(adjF,3),:));
        FN_prev = normalizerow(FN_prev);
        
        % reconnect faces
        F(F == vj) = vi;
        v = vi;
        tmp = [F, F(:,1)];
        [delF, ~] = find(diff(tmp,1,2) == 0); % find faces that to be deleted
        
        % compute dist(vi,optV) and dist(vj,optV)
        dist_v_vi = sqrt(sum((EOptPos(e,:) - V(vi,:)).^2));
        dist_v_vj = sqrt(sum((EOptPos(e,:) - V(vj,:)).^2));
        
        % move vertex i
        V(v,:) = EOptPos(e,:);
        
        
        % compute new normal
        for ii = 1:length(delF)
            FN_prev(adjF == delF(ii),:) = [];
            adjF(adjF == delF(ii)) = [];
        end
        FN_new = cross(V(F(adjF,1),:)-V(F(adjF,2),:), V(F(adjF,1),:) - V(F(adjF,3),:));
        FN_new = normalizerow(FN_new);
        
        % CHECK: adj face normals flip
        dotProd = sum(FN_new .* FN_prev,2);
        if any(dotProd < dotn_threshold) 
            fprintf('normal flipped\n');
            V = Vpre;
            F = Fpre;
            ECost(e) = inf;
            if stall > maxStall break; end
            continue;
        end
        
        % CHECK: triangle quality
        l0 = sqrt(sum((V(F(adjF,1),:)-V(F(adjF,2),:)).^2,2));
        l1 = sqrt(sum((V(F(adjF,2),:)-V(F(adjF,3),:)).^2,2));
        l2 = sqrt(sum((V(F(adjF,3),:)-V(F(adjF,1),:)).^2,2));
        x = (l0+l1+l2) ./ 2;
        delta = sqrt(x .* (x-l0) .* (x-l1) .* (x-l2));
        triQ = 4 * sqrt(3) * delta ./ (l0.^2 + l1.^2 + l2.^2);
        if any(triQ < triQ_threshold)
            fprintf('bad triangle quality\n');
            V = Vpre;
            F = Fpre;
            ECost(e) = inf;
            if stall > maxStall break; end
            continue;
        end
        
        % start post-collapsing
        numCollapse = numCollapse + 1;
        
        % adjacent vertex indices of v
        adjV = unique(F(adjF,:)); 
        adjV(adjV == v) = [];
        
        % delete face
        F(delF,:) = [];

        % reconnect edges
        assert(length(E) == length(ECost));
        assert(length(E) == size(EOptPos,1));
        assert(length(E) == length(ETag));
        E(E == vj) = vi;
        E(e,:) = [];
        EOptPos(e,:) = [];
        ECost(e) = [];
        ETag(e) = [];
        assert(length(E) == length(ECost));
        assert(length(E) == size(EOptPos,1));
        assert(length(E) == length(ETag));
        
        % update restriction
        R(vi,:) = (R(vi,:)*dist_v_vj + R(vj,:)*dist_v_vi) / (dist_v_vi + dist_v_vj);
        
        % update vertex quadric
        [ECost, KV, EOptPos] = updateEdgeCost(V,E,vi,vj,KV,ECost,EOptPos,method);
        
        stall = 0;
        
        %% check terminate
        if numCollapse == totalCollapse
            break;
        end
        if stall > 500
            break;
        end
    end
    
    [RV,IM,RVIdx] = remove_unreferenced(V,F);
    RF = IM(F);
    R = R(RVIdx,:);
end

function [E,ECost,EOptPos,ETag, KV] = edgeQuadric(V,F,bdWScale, method)

    E = edges(F);
    ECost = zeros(size(E,1),1);
     
    % get boundary edges/faces
    ET = edge_triangle_adjacency(F,E);
    [bdEIdx,~] = find(ET == -1);
    bdFIdx = ET(bdEIdx,1);
    
    % get edge tag
    ETag = zeros(size(E,1),1);
    ETag(bdEIdx) = 1; % 0: interior edge, 1: boundary edge
    
    % compute face area
    FN = faceNormals(V,F);
    FA = doublearea(V,F) / 2;
    FA(bdFIdx) = FA(bdFIdx) * bdWScale;

    % compute face quadrics
    KF = zeros(4,4,size(F,1));
    
    % for interior faces
    inFIdx = [1:size(F,1)];
    inFIdx(bdFIdx) = [];
    for ii = 1:length(inFIdx)
        fIdx = inFIdx(ii);
        d = - FN(fIdx,:) * V(F(fIdx,1),:)';
        p = [FN(fIdx,:), d]';
        KF(:,:,fIdx) = FA(fIdx) * (p * p'); % weight by FA
    end
    
    % for boundary faces
    for ii = 1:length(bdFIdx)
        eIdx = bdEIdx(ii);
        vi = E(eIdx,1);
        vj = E(eIdx,2);
        eVec = V(vi,:) - V(vj,:);
        eVec = eVec / norm(eVec);
        
        fIdx = bdFIdx(ii);
        fN = FN(fIdx,:);
        
        biN = cross(eVec, fN);
        biN = biN / norm(biN);
        
        d = - biN * V(vi,:)';
        p = [biN, d]';
        KF(:,:,fIdx) = FA(fIdx) * (p * p'); % weight by FA
    end
    
    % compute vertex quadrics
    adjFList = vertexFaceAdjacencyList(F);
    KV = zeros(4,4,size(V,1));
    for ii = 1:size(V,1)
        adjF = adjFList{ii};
        KV(:,:,ii) = sum(KF(:,:,adjF),3);
    end        
    
    % compute edge quadrics
    ECost = zeros(size(E,1),1);
    EOptPos = zeros(size(E,1),3);
    for e = 1:size(E,1)
        vi = E(e,1); 
        vj = E(e,2);
        Ke = KV(:,:,vi) + KV(:,:,vj); 
        
        [p, cost] = optPos(Ke, V(vi,:), V(vj,:));
        ECost(e) = cost;
        EOptPos(e,:) = p;
    end
    
    if strcmp(method, 'midPoint') 
        ECost = sqrt(sum((V(E(:,1),:) - V(E(:,2),:)).^2, 2));
    end
end

function [p, cost] = optPos(Ke, vip, vjp)
    % v'Kv = v'Av + 2b'v + c 
    A = Ke(1:3, 1:3);
    b = Ke(1:3, 4);
    c = Ke(4, 4);
    vK = A \ -b;
    vK = vK';

    costK = [vK,1] * Ke * [vK,1]';
    costi = [vip,1] * Ke * [vip,1]';
    costj = [vjp,1] * Ke * [vjp,1]';
    costMid = [(vip+vjp)./2,1]  * Ke * [(vip+vjp)./2,1]';

    [cost,idx] = min([costK, costi, costj, costMid]);
    posList = [vK; vip; vjp; (vip+vjp)./2];
    p = posList(idx,:);
end

function [ECost, K, EOptPos] = updateEdgeCost(V,E,vi,vj,K,ECost,EOptPos, method)
    K(:,:,vi) = K(:,:,vi) + K(:,:,vj);
    [adjE, ~] = find(E == vi);
    for ii = 1:length(adjE)
        eIdx = adjE(ii);
        Ke = K(:,:,E(eIdx,1)) + K(:,:,E(eIdx,2)); 
        [p, cost] = optPos(Ke, V(E(eIdx,1),:), V(E(eIdx,2),:));
        ECost(eIdx) = cost;
        EOptPos(eIdx,:) = p;
        if isnan(ECost(eIdx)) || isinf(ECost(eIdx))
            ECost(eIdx) = inf; % set inf or nan to inf
        end
    end
    K(:,:,vj) = nan;      
    
    if strcmp(method, 'midPoint')
        ECost(adjE) = sqrt(sum((V(E(adjE,1),:) - V(E(adjE,2),:)).^2,2));
    end
end