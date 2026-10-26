function [res] = RMS(X, Xhat)
%RMSE calculates RMS of X
    %X is array of scalars or vectors
arguments (Input)
    X
    Xhat
end

arguments (Output)
    res
end

errs = sqrt(sum((X-Xhat).^2,2));
res = sqrt(mean(errs.^2));

end