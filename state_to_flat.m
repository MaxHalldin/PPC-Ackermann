function y = state_to_flat(x,params)
%   Transforms the non-flat system to the flat one
%   x       -   state vector
%   params  -   constants
arguments (Input)
    x
    params
end

arguments (Output)
    y
end
    L = params{1};

    thetadot = x(5)/L*tan(x(4));

    y = zeros(6,1,'double');
    y(1) = x(1);
    y(2) = x(2);
    y(3) = x(5)*cos(x(3));
    y(4) = x(5)*sin(x(3));
    y(5) = x(6)*cos(x(3)) - x(5)*thetadot*sin(x(3));
    y(6) = x(6)*sin(x(3)) + x(5)*thetadot*sin(x(3));
end