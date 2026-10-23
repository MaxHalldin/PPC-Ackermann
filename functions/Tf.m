function res = Tf(x)
%APC scaling factor 
%   x in (-1,1)
%   x can be scalar or vector (or a matrix)
arguments (Input)
    x
end

arguments (Output)
    res
end
 
I = ones(size(x));
res = 0.5 * log( I+x ./ (I-x) );

end