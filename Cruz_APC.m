clear all; 
clc; 
close all;

%% 

A = 1;
c = 1; 
yr  = @(t)A*sin(c*t);
yp  = @(t)A*c*cos(c*t);
ref = {yr, yp};
tspan = [0; 7]; % Seconds

sampling_freq = 1e3; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

param = struct;

param.Phi = 0.2;
param.b0 = 1.85;
param.M = 0.12;
param.l = 0.28;

param.e_inf = 0.05;
param.lu = 1.5;
param.lambda = 2;
param.a = flip(poly(-1*param.lambda)); %Row vector % Need to have the negative here for the roots(lambdas).
param.rho_0 = 5;
param.rho_inf = param.lambda(1)*param.e_inf;
param.K = 1;
param.ubar = 2;
param.ebar = 2;

param.A = [0 1; 0 0];
param.B = [0; param.b0];
param.T = [1 0; param.Phi 1];
param.Tinv = [1 0; -param.Phi 1];

param.mu = 1;
param.L = 2.5;
param.k = [2*sqrt(3);6]; % values borrowed from article.
cond1 = ( 0<param.k(1) && param.k(1)<2*sqrt(param.L) ) && ...
        ( param.k(2) > param.k(1)^2/4 + 4*param.L/(param.k(2)^2) );
cond2 = ( param.k(1)>2*sqrt(param.L) && param.k(2) > 2*param.L );
assert(cond1||cond2);

ehat_0 = [-0.2; 0.3];
y_0 = [1.2; 1.5];
xer0 = [y_0; ehat_0; param.rho_0];

%%

Y = zeros(length(xer0),length(ts));
Y(:,1) = xer0;
U = zeros(2,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

euler = true;
tic
tot = length(ts);
if euler
    for i = 1:length(ts)-1
        % This is the forward Euler method. 
        [temp,U(:,i+1)] = URED_APC(ts(i),      Y(:,i),                  ref,param);
        Y(:,i+1) = Y(:,i) + dt*temp;
    end
else
    for i = 1:length(ts)-1
        % Runge-Kutta 4
        % Looks like it is about 3-4x slower than euler forward. Should
        % compare the score when i get it working later.
        [f1,U(i+1)] = URED_APC(ts(i),         Y(:,i),                 ref,param);
        [f2,~]      = URED_APC(ts(i)+1/3*dt,  Y(:,i)+dt/3*f1,         ref,param);
        [f3,~]      = URED_APC(ts(i)+2/3*dt,  Y(:,i)+dt*(-f1/3+f2),   ref,param);
        [f4,~]      = URED_APC(ts(i)+dt,      Y(:,i)+dt*(f1-f2+f3),   ref,param);
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
e = param.T*(Y(1:2,:)-[ref{1}(ts);ref{2}(ts)])-Y(3:4,:);
plot(ts,param.a*e)
plot(ts,Y(5,:),'r--')
plot(ts,-Y(5,:),'r--')
grid on
xlabel('t')
legend('True error', 's(e)', 'Error funnel')
title('error estimation')
hold off

figure
hold on
plot(ts,U(1,:),'r')
plot(ts,U(2,:))
plot(ts,param.ebar*ones(size(ts)),':r')
plot(ts,-param.ebar*ones(size(ts)),':r')
grid on
legend('u(t)','u_d(t)')
xlabel('t')
title('Control input')
hold off

figure
hold on
plot(ts,Y(3,:))
grid on
xlabel('t')
title('evolution of ehat_1')
hold off

figure
hold on
plot(ts,Y(4,:))
grid on
xlabel('t')
title('evolution of ehat_2')
hold off

%%
function [res,U] = URED_APC(t,y,ref,params)
    Y = y(1:2);
    ehat = y(3:4);
    rho = y(5);
   
    X = params.T*Y;
    X_ref = params.T*cellfun(@(c) c(t), ref).';
    e = X-X_ref;

    % URED
    % sigma(1) = e
    % sigma(2) = edot
    % z(1)(t) = ehat
    % z(2)(t) = ehatdot
    % Thus sigma(2) is the estimator error.
    sigma = ehat(1)-e(1);
    
    varphi = [ abs(sigma)^(1/2)*sign(sigma) + params.mu*abs(sigma)^(3/2)*sign(sigma);
            0.5*sign(sigma) + 2*params.mu*sigma + 1.5*params.mu^2*abs(sigma)^2*sign(sigma) ];
    
    ehatdot = params.A*ehat-params.k.*varphi;
    
    % Feedback control
    % 1. Linear filter
    % s = params.a*sat([params.T*e(1);ehat(2)], params.ebar);
    s = params.a*sat([e(1); ehatdot(1)], params.ebar);

    % 2. APC
    xi = s/rho;
    ud = -params.K*Jf(xi).*Tf(xi);
    u = sat(ud, params.ubar);
    
    rhodot = (u-ud)/xi - params.lu*(rho - params.rho_inf);

    phi = [-params.Phi*Y(1); -params.b0*( params.M*Y(1)^3 + params.l*Y(1) )];
    xdot = params.A*X + phi + params.B*u;
    ydot = params.Tinv*xdot;


    res = [ydot; ehatdot; rhodot];
    U = [u; ud];
end

%% Helper functions

function res = Tf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    % assert(isbetween(x,-1,1), "Argument is out of bounds.");
    I = ones(size(x));
    res = 0.5 * log( I+x ./ (I-x) );
end

function res = Jf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    % assert(isbetween(x,-1,1), "Argument is out of bounds.");
    I = ones(size(x));
    res = 1./(I-x.^2);
end

function res = sat(x, xbar)
    % x can be vector or scalar. If x is a scalar, xbar must also be. If x 
    % is a vector, then xbar can either be a scalar or a vector of same 
    % size as x.
    res = sign(x).*min(abs(x), xbar);
end