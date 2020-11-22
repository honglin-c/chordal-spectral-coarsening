function matrixChecks(X, M, matStr)
    if(~exist('matStr','var'))
        matStr = 'matrix';
    end

    if max(max(abs(X - X'))) > 1e-5
        warning(strcat(matStr, " is not symmetric"))
    end
    if max(abs(X*ones(size(X,1),1))) > 1e-5
        warning(strcat(matStr, " violates null-space constraint"))
    end
    if size(X,1) < 20
      [eVal, ~] = eigsReal(X,M,size(X,1));
    else
      [eVal, ~] = eigsReal(X,M,20);
    end
    
    if min(eVal) < -1e-5
        warning(strcat(matStr, " is not PSD, min eVal = ", num2str(min(eVal))))
    end
end
        