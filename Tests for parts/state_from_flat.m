function x = state_from_flat(y,L)
%   Transforms the flat system to non-flat
%   x       -   state vector
%               -> "x1" = "y1"
%               -> "x2" = "y2"
%               -> "x3" = "theta"
%               -> "x4" = "phi"
%               -> "x5" = "u1"
%               -> "x6" = "w"
%   params  -   constants
arguments (Input)
    y  
    L
end

arguments (Output)
    x
end
    x = zeros(6,1,'double');
    x(1) = y(1);
    x(2) = y(2);
    x(3) = atan(y(4)/y(3));
    x(5) = sqrt( y(3)^2 + y(4)^2 );
    x(4) = atan( L*(y(3)*y(6) - y(5)*y(4)) / x(5)^3 );
    x(6) = (y(3)*y(5) + y(6)*y(4))/x(5);
end