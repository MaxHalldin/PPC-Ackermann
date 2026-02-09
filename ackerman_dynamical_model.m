clear all; clc; close all;

y1t  = @(t)sin(t);
y2t  = @(t)cos(t);
y1tp = @(t)cos(t);
y2tp = @(t)-sin(t);
y1tb = @(t)-sin(t);
y2tb = @(t)-cos(t);
y1tt = @(t)-cos(t);
y2tt = @(t)sin(t);
traj = {y1t y2t y1tp y2tp y1tb y2tb y1tt y2tt};

function xdot = ODE(t, x, params, ref)
    % Function based on the psuedocode of Arcuri, et al.
    [x1, x2, theta, phi, u1, w] = deal(x(1),x(2),x(3),x(4),x(5),x(6));
    [L, k0, k1, k2] = deal(params(1),params(2),params(3),params(4));
    [x1ref, x2ref, x1pref, x2pref, x1bref, x2bref, x1tref, x2tref] = ref{:};

    x1dot = u1*cos(theta);
    x2dot = u1*sin(theta);
    thetadot = u1/L*tan(phi);
    u1dot = w;

    x1ddot = u1dot*cos(theta) - u1*thetadot*sin(theta);
    x2ddot = u1dot*sin(theta) + u1*thetadot*cos(theta);

    e1 = x1-x1ref(t);          e2 = x2-x2ref(t);
    e1p = x1dot-x1pref(t);     e2p = x2dot-x2pref(t);
    e1b = x1ddot-x1bref(t);    e2b = x2ddot-x2bref(t);

    a1 = x1tref(t) - k2*e1b - k1*e1p - k0*e1;
    a2 = x2tref(t) - k2*e2b - k1*e2p - k0*e2;

    C1 = -3*w*u1/L*tan(phi)*sin(theta) - u1^3/L*tan(phi)^2*cos(theta);
    C2 = 3*w*u1/L*tan(phi)*cos(theta) - u1^3/L*tan(phi)^2*sin(theta);

    k = u1^2/L*sec(phi)^2;
    k = max([k 1e-3]);
    
    v1 = (a1-C1)*cos(theta) + (a2-C2)*sin(theta);
    v2 = ((a2-C2)*cos(theta) - (a1-C1)*sin(theta))/k ;

    wdot = v1;
    phidot = v2;
    xdot = [x1dot; x2dot; thetadot; phidot; u1dot; wdot];
end

function testTrack(ref)
    tspan = [0 20];
    params = [1.0 3.0 6.0 4.5];
    y0 = [1 1 0.2 0.05 0.15 0];
    
    [t, X] = ode45(@(t,y)ODE(t,y,params,ref), tspan, y0);
    
    e1 = X(:,1)-ref{1}(t);
    e2 = X(:,2)-ref{2}(t);
    
    f1 = figure;
    plot(X(:,1), X(:,2),ref{1}(t),ref{2}(t),':')
    grid on
    xlabel('x1')
    ylabel('x2')
    legend('estimated path', 'true path')
    title('Path tracking')
    f2 = figure;
    plot(t,e1)
    hold on
    plot(t,e2)
    grid on
    xlabel('time')
    ylabel('error')
    legend('e_1', 'e_2')
    title('Tracking error')
    hold off
end

testTrack(traj)