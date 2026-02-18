function [dydt, inputs] = ydotPPC(t, Y, ref, params, rhos)
%   ydot outputs the first time derivative of the flat Ackermann system state
%   vector. Also outputs the motor force and rate of change of steering angle.
%   t       - 1x1 time. Used for trajectory functions.
%   Y       - 6x1 State vector of flat system. 
%   ref     - 6x1 cell array containing the reference trajectory functions.
%   params  - 7x1 cell array of vehicle parameters. (length & mass of axles)
%   rhos    - 3x1 cell array containing performance functions.
arguments (Input)
    t
    Y
    ref
    params
    rhos
end

arguments (Output)
    dydt
    inputs
end
    [L, M1, M2, k] = params{:};
    [rho1, rho2, rho3] = rhos{:};
    R = cellfun(@(c) c(t), ref);
    
    % maybe a1&2 should be in dydt.
    a1  = -k(1)*Tf( (Y(1:2) - R(1:2).') / rho1(t) );
    a2  = -k(2)*Tf( (Y(3:4) - a1)       / rho2(t) );
    u   = -k(3)*Tf( (Y(5:6) - a2)       / rho3(t) );

    A = [zeros(4,2) eye(4,4); zeros(2,6)];
    B = [zeros(4,2); eye(2,2)];

    dydt = A*Y+B*u;

    %%
    u1 = sqrt(Y(3)^2+Y(4)^2);
    phi = atan((Y(3)*Y(6)-Y(4)*Y(5))/(u1^3)*L);
    u1dot = (Y(3)*Y(5)+Y(4)*Y(6))/u1;

    phidot = (L*sqrt(Y(3)^2 + Y(4)^2) * (u(2) * Y(3) - Y(4) * u(1)) *...
             (Y(3)^2 + Y(4)^2) - 3*(Y(5) * Y(3) + Y(6) * Y(4)) *...
             (Y(6) * Y(3) - Y(5) * Y(4))) / ((Y(3)^2 + Y(4)^2)^3 + L^2 *...
             (Y(6) * Y(3) - Y(5)*Y(4))^2);
    F = (M1+M2+M2*tan(phi)^2)*u1dot + M2*phidot*sin(phi)*u1/(cos(phi)^3);

    inputs = [phidot; F];
    % disp("hit")
end

function res = Tf(x)
    assert((x(1) > -1 && x(1) < 1), "Argument 1 %d is out of bounds", x(1));
    assert((x(2) > -1 && x(2) < 1), "Argument 2 %d is out of bounds", x(2));
    res = arrayfun(@(z) log((1+z)/(1-z)), x );
end