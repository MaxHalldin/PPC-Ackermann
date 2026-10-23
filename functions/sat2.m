function res = sat2(x,xbar)
%sat Non-symmetrical saturation function
%   x       - scalar
%   xbar    - vector of upper and lower saturation of x
%   xbar(1)<xbar(2)
arguments (Input)
    x
    xbar
end

arguments (Output)
    res
end

    res = x + (xbar(1)-x)*(xbar(1)<x) + (xbar(2)-x)*(x<xbar(2));
end