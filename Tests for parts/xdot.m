function [dydt, inputs] = xdot(t, X, ref, params)
%   xdot outputs the first time derivative of the Ackermann system state
%   vector. Also outputs the motor force and rate of change of steering angle.
%   t       - 1x1 time. Used for trajectory functions.
%   X       - 6x1 State vector of system. 
%   ref     - 6x1 cell array containing the reference trajectory.
%   params  - 3x1 cell array of vehicle parameters. (length & mass of axles)
arguments (Input)
    t 
    X
    ref
    params
end

arguments (Output)
    dydt
    inputs
end

    [x1, x2, theta, phi, u1, w] = deal(X(1),X(2),X(3),X(4),X(5),X(6));
    [L, M1, M2] = deal(params{:});
    [y1ref, y2ref, y1pref, y2pref, y1bref, y2bref] = ref{:};

    % Calculate states needed to transition to flat system
    % TODO: add upper and lower bounds for these parameters.
    thetadot = u1*tan(phi)/L;
    u1dot = w;

    % Map system state to flat system
    y1 = x1; 
    y2 = x2;
    y1dot = u1*cos(theta);
    y2dot = u1*sin(theta);
    y1ddot = w*cos(theta) - u1*thetadot*sin(theta);
    y2ddot = w*sin(theta) + u1*thetadot*cos(theta);

    % Express flat system in matrix form
    Y = [y1 y2 y1dot y2dot y1ddot y2ddot].';
    R = [y1ref(t) y2ref(t) y1pref(t) y2pref(t) y1bref(t) y2bref(t)].';
    A = [zeros(4,2) eye(4,4); zeros(2,6)];
    B = [zeros(4,2); eye(2,2)];

    % Poles placement based on Franch, et al.
    p = [-200 -200 -100 -100 -2 -2]; 
    K = place(A,B,p);
    V = K*(R-Y);

    y1dddot = V(1);
    y2dddot = V(2);

    phidot = (L*sqrt(y1dot^2 + y2dot^2) * ...
             (y2dddot * y1dot - y2dot * y1dddot) *...
             (y1dot^2 + y2dot^2) - 3*...
             (y1ddot * y1dot + y2ddot * y2dot) *...
             (y2ddot * y1dot - y1ddot * y2dot)) /...
             ((y1dot^2 + y2dot^2)^3 + L^2 *...
             (y2ddot * y1dot - y1ddot*y2dot)^2);
    F = (M1 + M2 + M2*tan(phi)^2)*u1dot +... 
        M2*phidot*sin(phi)*u1/(cos(phi)^3);
    
    a = -u1dot*u1/L*tan(phi)*sin(theta) - u1^3/L*tan(phi)^2*cos(theta);
    b = u1dot*u1/L*tan(phi)*cos(theta) - u1^3/L*tan(phi)^2*sin(theta);
    wdot = (y1dddot-a)*cos(theta) + (y2dddot-b)*sin(theta);

    inputs = [phidot; F];
    dydt = [y1dot; y2dot; thetadot; phidot; u1dot; wdot];
end