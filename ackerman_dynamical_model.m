clear all; clc; close all;

% Define trajectory. For now circle.
y1t  = @(t)sin(t);
y2t  = @(t)cos(t);
y1tp = @(t)cos(t);
y2tp = @(t)-sin(t);
y1tb = @(t)-sin(t);
y2tb = @(t)-cos(t);
traj = {y1t y2t y1tp y2tp y1tb y2tb};

function [dydt, inputs] = ydot(t, x, params, ref)
    [x1, x2, theta, phi, u1, w] = deal(x(1),x(2),x(3),x(4),x(5),x(6));
    [L, M1, M2] = deal(params(1));
    [y1ref, y2ref, y1pref, y2pref, y1bref, y2bref] = ref{:};

    % Calculate states needed to transition to flat system
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

ref = traj;
tspan = [0 20];
params = [1.0 2.0 2.0];
y0 = [1 1 0.2 0.05 0.15 0];

[t, Y] = ode45(@(t,y)ydot(t,y,params,ref), tspan, y0);

figure
hold on
plot(Y(:,1), Y(:,2))
plot(ref{1}(t),ref{2}(t),':')
grid on
xlabel('x1')
ylabel('x2')
legend('Estimated path', 'True path')
title('Path tracking')
hold off

figure
hold on
for i=1:6
    e = Y(:,i)-ref{i}(t);
    plot(t,e)
end
xlabel('Time')
ylabel('Error')
legend('x', 'y', 'xdot', 'ydot', 'xddot', 'yddot')
title('Position error error')
grid on
hold off

figure
inputs = zeros(length(t),2);
for k = 1:length(inputs)
    [~,inputs(k,:)] = ydot(t(k),Y(k,:),params,ref);
end
% Using log plot due to early values are quite large, though also include
% zeros
% TODO: add upper and lower bounds for F and Phidot.
semilogy(t,inputs(:,1),t,inputs(:,2)) 
title('Position error error')
legend("phidot", "F")
grid on