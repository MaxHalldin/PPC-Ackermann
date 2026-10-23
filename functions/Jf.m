function [res] = Jf(x)
%Jf jacobian of APC scalar function
%   x in (-1,1)
%   x can be scalar or vector (or a matrix)

arguments (Input)
    x
end

arguments (Output)
    res
end

    I = ones(size(x));
    res = 1./(I-x.^2);
end