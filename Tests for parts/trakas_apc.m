clear all; 
clc; 
close all;

%% I suspect that my problem is that I need to sample faster than ode45 
% does for integrating

A = 1;
c = 1; 
yr  = @(t)A*sin(c*t);
yp  = @(t)A*c*cos(c*t);
ref = {yr, yp};
tspan = [0; 50]; % Seconds

sampling_freq = 1e4; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

param.Phi = 0.2;
param.b0 = 1.85;
param.M = 0.12;
param.L = 0.28;

param.e_inf = 0.05;
param.ebar = 2;
param.lambda = 2;
param.a = flip(poly(-1*param.lambda)); %Row vector
param.gamma = [1.1; 1.0];


param.lu = 1.5;
param.rho_0 = 5;
param.rho_inf = prod(param.lambda)*param.e_inf;
param.k = 1;
param.ubar = 2;

param.A = [0 1; 0 0];
param.B = [0; param.b0];
param.T = [1 0; param.Phi 1];
param.Tinv = [1 0; -param.Phi 1];

%%
% If the initial guess is correct than it works well, which suggests to me
% that there is some set for which the PPC observer works correctly.
y_0 = [1.2; 1.5];
x0 = [y_0; param.rho_0];

Y = zeros(length(x0),length(ts));
Y(:,1) = x0;
U = zeros(1,length(ts));

euler = true;
tic
tot = length(ts);
if euler
    for i = 1:length(ts)-1
        % This is the forward Euler method. 
        [f,U(i+1)] = ydotPPCtest(ts(i),Y(:,i),ref,param);
        Y(:,i+1) = Y(:,i) + dt*f;
    end
else
    for i = 1:length(ts)-1
        % Runge-Kutta 4
        % Looks like it is about 3-4x slower than euler forward. Should
        % compare the score when i get it working later.
        % This should only really be used for simulation, as we want a
        % fixed sample rate for the robot implementation.
        [f1,U(i+1)] = ydotPPCtest(ts(i),         Y(:,i),                 ref,param);
        [f2,~]      = ydotPPCtest(ts(i)+1/3*dt,  Y(:,i)+dt/3*f1,         ref,param);
        [f3,~]      = ydotPPCtest(ts(i)+2/3*dt,  Y(:,i)+dt*(-f1/3+f2),   ref,param);
        [f4,~]      = ydotPPCtest(ts(i)+dt,      Y(:,i)+dt*(f1-f2+f3),   ref,param);
        Y(:,i+1)    = Y(:,i) + dt/8*(f1+3*f2+3*f3+f4);
    end
end
toc

%% PLOTTING

figure
hold on
plot(ts,Y(1,:))
plot(ts,ref{1}(ts),':')
grid on
xlabel('t')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nMSE = %0.3e', mean( (Y(1,:)-ref{1}(ts)).^2 )))
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts),':')
grid on
xlabel('t')
legend('Estimated velocity', 'True velocity')
title(sprintf('Velocity tracking\nMSE = %0.3e', mean( (Y(2,:)-ref{2}(ts)).^2 )))
hold off

figure
hold on
plot(ts,Y(1,:)-ref{1}(ts))
e = param.T*(Y(1:2,:)-[ref{1}(ts);ref{2}(ts)]);
plot(ts,param.a*e)
plot(ts,Y(3,:),'r--')
plot(ts,-Y(3,:),'r--')
grid on
xlabel('t')
legend('True error', 's(e)', 'Error funnel')
title('error estimation')
hold off

figure
hold on
plot(ts,U,'r')
plot(ts,param.ebar*ones(size(ts)),':r')
plot(ts,-param.ebar*ones(size(ts)),':r')
grid on
ylabel('u(t)')
xlabel('t')
title('Control input')
hold off

%%
function [res,u] = ydotPPCtest(t,y,ref,params)
    Y = y(1:2);
    rho = y(3);

    X = params.T*Y;

    % Need to transform coordinated to state variables. (Liu, et. al.)
    Y_ref = cellfun(@(c) c(t), ref).';

    X_ref = params.T*Y_ref;

    %This phi is in accordance with Trakas, et. al. but not consistant with
    %the source thas is cited.
    phi = [-params.Phi*Y(1); -params.b0*( params.M*Y(2) + params.L*Y(1) )];

    e = X-X_ref;

    % Feedback control
    % 1. Linear filter
    s = params.a*sat(e, params.ebar);

    % 2. APC
    xi = s/rho;
    ud = -params.k*Jf(xi).*Tf(xi);
    u = sat(ud, params.ubar);

    rhodot = -params.lu*(rho - params.rho_inf) + (u-ud)/xi;
    xdot = params.A*X + phi + params.B*u;
    ydot = params.Tinv*xdot;

    res = [ydot; rhodot];
end

%% Helper functions

function res = Tf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    I = ones(size(x));
    res = 0.5 * log( I + x ./ (I - x) );
end

function res = Jf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    I = ones(size(x));
    res = 1./(I-x.^2);
end

function res = sat(x, xbar)
    % x can be vector or scalar. If x is a scalar, xbar must also be. If x 
    % is a vector, then xbar can either be a scalar or a vector of same 
    % size as x.
    res = sign(x).*min(abs(x), xbar);
end