clear all; 
clc; 
close all;

% Define trajectory. For now circle.
c = pi/10;
y1t  = @(t)sin(c*t);
y2t  = @(t)cos(c*t);
y1tp = @(t)c*cos(c*t);
y2tp = @(t)-c*sin(c*t);
y1tb = @(t)-c^2*sin(c*t);
y2tb = @(t)-c^2*cos(c*t);
traj = {y1t y2t y1tp y2tp y1tb y2tb};

function y = state_to_flat(x,params)
    [x1, x2, theta, phi, u1, w] = deal(x(1),x(2),x(3),x(4),x(5),x(6));
    [L, ~, ~] = deal(params{:});

    thetadot = u1/L*tan(phi);

    y1 = x1;
    y2 = x2;
    y1dot = u1*cos(theta);
    y2dot = u1*sin(theta);
    y1ddot = w*cos(theta) - u1*thetadot*sin(theta);
    y2ddot = w*sin(theta) + u1*thetadot*sin(theta);

    y = [y1; y2; y1dot; y2dot; y1ddot; y2ddot];
end 

ref = traj;
tspan = [0; 20];
params = {1.0 2.0 2.0};
x0 = [1; 1; 0.2; 0.05; 0.15; 0];

tic
[t, Y] = ode113(@(t,y)xdot(t,y,ref,params), tspan, x0);
% [t, Y] = ode15s(@(t,y)ydot(t,y,ref,params), tspan, y0);
toc

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
e = zeros(size(Y));
for i=1:2
    e(:,i) = Y(:,i)-ref{i}(t);
end
plot(t,e(:,1:2))
xlabel('Time')
ylabel('Error')
legend('x', 'y')
title('Position error error')
grid on

figure
inputs = zeros(length(t),2);
for k = 1:length(inputs)
    [~,inputs(k,:)] = ydot(t(k),Y(k,:).',ref,params);
end
plot(t,inputs(:,1),t,inputs(:,2)) 
title('Position error error')
legend("phidot", "F")
grid on