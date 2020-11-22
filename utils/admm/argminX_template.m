function X = argminX_template(z,u,rho,dLHS,data)
    RHS = [data.Ea; data.P.z2xc*(z- u); sparse(data.nG,1)];
    
    % solve KKT
    Xmu = dLHS \ RHS;

    % get y
    x = Xmu(1:data.nx);
    mu1 = Xmu(data.nx+1 : data.nx+data.nxc);
    y = z - u + 1/rho * data.P.z2xc' * mu1;
    
    assert( all(abs(data.P.x2xc*x - data.P.z2xc*y) < 1e-7)) % check constraint
    
    X = [x;y];
end