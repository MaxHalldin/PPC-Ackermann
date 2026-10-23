function res = sat(x,xbar)
%sat Saturation function
%   x can be vector or scalar. If x is a scalar, xbar must also be. If x 
% is a vector, then xbar can either be a scalar or a vector of same size 
% as x.
arguments (Input)
    x
    xbar
end

arguments (Output)
    res
end

    res = sign(x).*min(abs(x), xbar);
end