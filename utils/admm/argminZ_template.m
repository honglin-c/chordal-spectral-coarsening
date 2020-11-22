function z = argminZ_template(xy,u,rho,data)
    y = xy(data.nx+1:end);
    Y = data.P.z2Z * y;
    U = data.P.z2Z * u;
    
    z = zeros(data.nz,1);
    zIdx = 0;
    for k = 1:length(data.cliques)
        ncl = length(data.cliques{k});
        idx = data.vecCliques{k};
        toproject = reshape(Y(idx) + U(idx), ncl, ncl);
        mass = data.Mc(data.cliques{k}, data.cliques{k});
        Zk = projectPSD(toproject, mass, 0); % project sub-matrix to PSD 
        zk = data.P.Zk2zk{k} * Zk(:);
        z(zIdx+1 : zIdx+length(zk)) = zk;
        zIdx = zIdx + length(zk);
    end
    
end