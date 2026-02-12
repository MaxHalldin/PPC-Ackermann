function [dydt, inputs] = ydot(t, Y, ref, params)
%   ydot outputs the first time derivative of the flat Ackermann system state
%   vector. Also outputs the motor force and rate of change of steering angle.
%   t       - 1x1 time. Used for trajectory functions.
%   Y       - 6x1 State vector of flat system. 
%   ref     - 6x1 cell array containing the reference trajectory.
%   params  - 3x1 cell array of vehicle parameters. (length & mass of axles)
arguments (Input)
    t
    Y
    ref
    params
end

arguments (Output)
    dydt
    inputs
end

    [y1dot, y2dot, y1ddot, y2ddot] = deal(Y(3),Y(4),Y(5),Y(6));
    [y1ref, y2ref, y1pref, y2pref, y1bref, y2bref] = ref{:};
    [L, M1, M2] = deal(params{:});

    R = [y1ref(t) y2ref(t) y1pref(t) y2pref(t) y1bref(t) y2bref(t)].';
    A = [zeros(4,2) eye(4,4); zeros(2,6)];
    B = [zeros(4,2); eye(2,2)];

    % Poles placement based on Franch, et al.
    p = [-200 -200 -100 -100 -2 -2]; 
    K = place(A,B,p);
    V = K*(R-Y);

    y1dddot = V(1);
    y2dddot = V(2);

    u1 = sqrt(y1dot^2+y2dot^2);
    phi = atan((y1dot*y2ddot-y2dot*y1ddot)/(u1^3)*L);
    u1dot = (y1dot*y1ddot+y2dot*y2ddot)/u1;

    phidot = (L*sqrt(y1dot^2 + y2dot^2) * ...
             (y2dddot * y1dot - y2dot * y1dddot) *...
             (y1dot^2 + y2dot^2) - 3*...
             (y1ddot * y1dot + y2ddot * y2dot) *...
             (y2ddot * y1dot - y1ddot * y2dot)) /...
             ((y1dot^2 + y2dot^2)^3 + L^2 *...
             (y2ddot * y1dot - y1ddot*y2dot)^2);
    F = (M1 + M2 + M2*tan(phi)^2)*u1dot +... 
        M2*phidot*sin(phi)*u1/(cos(phi)^3);

    dydt = [Y(3:end); V];
    inputs = [phidot; F];
end