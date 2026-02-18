function x = state_from_flat(y,params)
%   Transforms the flat system to non-flat
%   x       -   state vector
%   params  -   constants
arguments (Input)
    y  
    params
end

arguments (Output)
    x
end
    L = params{1};

    x = zeros(6,1,'double');
    x(1) = y(1);
    x(2) = y(2);
    x(3) = atan(y(4)/y(3));
    x(5) = sqrt( y(3)^2 + y(4)^2 );
    x(4) = atan( L*(y(3)*y(6) - y(5)*y(4)) / x(5)^3 );
    x(6) = (y(3)*y(5) + y(6)*y(4))/x(5);
end