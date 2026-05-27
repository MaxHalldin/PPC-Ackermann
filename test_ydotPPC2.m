clear all; 
clc; 
close all;

% Define trajectory. For now circle.
A = 10;
c = pi/10;
y1t  = @(t)A*sin(c*t);
y2t  = @(t)A*cos(c*t);
y1tp = @(t)A*c*cos(c*t);
y2tp = @(t)A*-c*sin(c*t);
y1tb = @(t)A*-c^2*sin(c*t);
y2tb = @(t)A*-c^2*cos(c*t);
ref = {y1t y2t y1tp y2tp y1tb y2tb};

tspan = [0; pi/(6*c)];
dt = 1/100;
L = 1;
M1 = 1;
M2 = 1;
k = 1;
lambda = [-1 -2 ].';
ubar = 2e2*[1; 1];
gammas = [-1 -1.1 -1.3].';

params = {L, M1, M2, k, lambda, ubar, gammas};
x0 = [ref{1}(0)+0.01; ref{2}(0)+0.01; 2; 0; 0; 0];
y0 = state_to_flat(x0,params);

ehat0 = [0.5 0.5 0.5 0.5 0.5 0.5].';
einf = 1e-1;
ebar = 2*[1; 1; 1];

lu = [1.0; 1.0];
lo = [1.0; 1.0];
r_inf = einf*lambda(1);
r0 = 100*[1; 1];
rho0 = 100*[1; 1];
rho_inf = 1*[1; 1];

eparams = {ebar};
rhoparams = {r0, r_inf, rho_inf, lu, lo};

xer0 = [x0; ehat0; rho0];

tic
[t, Y] = ode45(@(t,y)ydotPPC2(t,y,ref,params,rhoparams,eparams), tspan, xer0);
% t = tspan(1):dt:tspan(2);
% 
% dydt = zeros(14,length(t));
% Y = zeros(14,length(t));
% Y(:,1) = xer0;
% 
% inputs = zeros(2,length(t));
% 
% for i = 1:length(t)
%     [dydt(:,i), inputs(:,i)] = ydotPPC2(t(i),Y(:,i),ref, params, rhoparams, eparams);
%     Y(:,i+1) = Y(:,i)+dydt(:,i)*dt;
% end
toc

%% PLOTTING
figure
hold on
plot(ref{1}(t),ref{2}(t),':')
plot(Y(:,1), Y(:,2))
grid on
xlabel('x1')
ylabel('x2')
legend('Estimated path', 'True path')
title('Path tracking')
hold off

figure
hold on
% plot(ref{1}(t),ref{2}(t),':')
plot(Y(:,13))
plot(Y(:,14))
grid on
ylabel('rho')
legend('rho1', 'rho1')
title('performance function')
hold off

% e = zeros(size(Y(:,1:6)));
% flat_state = zeros(size(Y(:,1:6)));
% for i=1:length(t)
%     flat_state(i,:) = state_to_flat(Y(i,1:6),params);
% end
% 
% for i=1:6
%     e(:,i) = flat_state(:,i)-ref{i}(t);
% end
% figure
% title("Tracking errors")
% ax1=subplot(3,1,1);
% hold on
% plot(t,e(:,1),t,e(:,2));
% plot(t,rho1(t), 'c--')
% plot(t,-rho1(t), 'c--')
% xticks([])
% legend('x', 'y')
% hold off
% ax2=subplot(3,1,2);
% hold on
% plot(t,e(:,3),t,e(:,4));
% plot(t,rho2(t), 'c--')
% plot(t,-rho2(t), 'c--')
% xticks([])
% ylabel('Error')
% legend('vx', 'vy')
% hold off
% ax3=subplot(3,1,3);
% hold on
% plot(t,e(:,5),t,e(:,6));
% plot(t,rho3(t), 'c--')
% plot(t,-rho3(t), 'c--')
% legend('ax', 'ay')
% hold off
% ax2.Position(2)=ax3.Position(2)+ax3.Position(4);
% ax1.Position(2)=ax2.Position(2)+ax2.Position(4);
% linkaxes([ax1 ax2 ax3],'x')


figure
inputs = zeros(length(t),2);
for k = 1:length(inputs)
    [~,inputs(k,:)] = ydotPPC2(t(k),Y(k,:).',ref,params,rhoparams,eparams);
end
% Using log plot due to early values are quite large, though also include
% zeros
% semilogy(t,inputs(:,1),t,inputs(:,2)) 
plot(t,inputs(:,1),t,inputs(:,2)) 
title('Position error')
legend("phidot", "F")
grid on