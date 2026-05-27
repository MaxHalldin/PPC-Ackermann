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
ref = {y1t y2t y1tp y2tp y1tb y2tb};

tspan = [0; 20];
params = {1.0, 2.0, 2.0, [40,30,100]};
x0 = [0; 1; 0.1; 0; 0; 0];
y0 = state_to_flat(x0,params);

l = 1;
rho_inf = 0.2;
rho0 = 50;

rho1 = @(t)(rho0*1-rho_inf)*exp(-1*t) + rho_inf;
rho2 = @(t)(rho0*3-rho_inf)*exp(-2*t) + rho_inf;
rho3 = @(t)(rho0/50-rho_inf)*exp(-0.5*l*t) + rho_inf;
rhos = {rho1, rho2, rho3};

tic
[t, Y] = ode45(@(t,y)ydotPPC(t,y,ref,params,rhos), tspan, x0);
toc

%% PLOTTING
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


e = zeros(size(Y));
flat_state = zeros(size(Y));
for i=1:length(t)
    flat_state(i,:) = state_to_flat(Y(i,:),params);
end

for i=1:6
    e(:,i) = flat_state(:,i)-ref{i}(t);
end
figure
title("Tracking errors")
ax1=subplot(3,1,1);
hold on
plot(t,e(:,1),t,e(:,2));
plot(t,rho1(t), 'c--')
plot(t,-rho1(t), 'c--')
xticks([])
legend('x', 'y')
hold off
ax2=subplot(3,1,2);
hold on
plot(t,e(:,3),t,e(:,4));
plot(t,rho2(t), 'c--')
plot(t,-rho2(t), 'c--')
xticks([])
ylabel('Error')
legend('vx', 'vy')
hold off
ax3=subplot(3,1,3);
hold on
plot(t,e(:,5),t,e(:,6));
plot(t,rho3(t), 'c--')
plot(t,-rho3(t), 'c--')
legend('ax', 'ay')
hold off
ax2.Position(2)=ax3.Position(2)+ax3.Position(4);
ax1.Position(2)=ax2.Position(2)+ax2.Position(4);
linkaxes([ax1 ax2 ax3],'x')


figure
inputs = zeros(length(t),2);
for k = 1:length(inputs)
    [~,inputs(k,:)] = ydotPPC(t(k),Y(k,:).',ref,params,rhos);
end
% Using log plot due to early values are quite large, though also include
% zeros
% semilogy(t,inputs(:,1),t,inputs(:,2)) 
plot(t,inputs(:,1),t,inputs(:,2)) 
title('Position error')
legend("phidot", "F")
grid on