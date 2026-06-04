function y = state_to_flat(x,L)
%   Transforms the non-flat system to the flat one
%   x - state vector
%    -> "x1" = "y1"
%    -> "x2" = "y2"
%    -> "x3" = "theta"
%    -> "x4" = "phi"
%    -> "x5" = "u1"
%    -> "x6" = "w" or "u1dot"
%   params  -   constants
arguments (Input)
    x
    L
end

arguments (Output)
    y
end

    thetadot = x(5)/L*tan(x(4));

    y = zeros(6,1,'double');
    y(1) = x(1);
    y(2) = x(2);
    y(3) = x(5)*cos(x(3));
    y(4) = x(5)*sin(x(3));
    y(5) = x(6)*cos(x(3)) - x(5)*thetadot*sin(x(3));
    y(6) = x(6)*sin(x(3)) + x(5)*thetadot*sin(x(3));
end