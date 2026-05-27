clear all; 
clc; 
close all;

%% 

A = 1;
c = pi/10;
y1t  = @(t)A*sin(c*t);
y2t  = @(t)A*cos(c*t);
y1tp = @(t)A*c*cos(c*t);
y2tp = @(t)A*-c*sin(c*t);
y1tb = @(t)A*-c^2*sin(c*t);
y2tb = @(t)A*-c^2*cos(c*t);
ref = {y1t y2t y1tp y2tp y1tb y2tb};

tspan = [0; 20]; % Seconds
sampling_freq = 1e3; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

param = struct;

param.l = 0.5;
param.M1 = 2;
param.M2 = 2;
param.Fbar = 10; % Placeholder
param.phidotbar = 10; % Placeholder

param.e_inf = 0.05*[1 1].';
param.lu = 1.5;
param.lambda = [2 1].';
param.a = flip(poly(-1*param.lambda)); %Row vector % Need to have the negative here for the roots(lambdas).
param.rho_0 = 5*[1 1].';
param.rho_inf = param.lambda(1)*param.e_inf;
param.K = 1.2;
param.ubar = 5;
param.ebar = 5;

param.A = [zeros(4,2) eye(4,4); zeros(2,6)];
param.B = [zeros(4,2); eye(2,2)];

param.mu = 1;
param.L = 2.5;
param.k = [2*sqrt(3);6]; % values borrowed from article.
cond1 = ( 0<param.k(1) && param.k(1)<2*sqrt(param.L) ) && ...
        ( param.k(2) > param.k(1)^2/4 + 4*param.L/(param.k(2)^2) );
cond2 = ( param.k(1)>2*sqrt(param.L) && param.k(2) > 2*param.L );
assert(cond1||cond2);

y0 = [-0.01; 0.99; 0.1; 0.01; 0; 0];

ehat_0 = [0 0 0 0 0 0].';
xer0 = [y0; ehat_0; param.rho_0];

%%

Y = zeros(length(xer0),length(ts));
Y(:,1) = xer0;
U = zeros(4,length(ts));
inputs = zeros(2,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

euler = true;
tic
tot = length(ts);
if euler
    for i = 1:length(ts)-1
        % This is the forward Euler method. 
        [temp,U(:,i+1),inputs(:,i+1)] = URED_APC(ts(i),      Y(:,i),                  ref,param);
        Y(:,i+1) = Y(:,i) + dt*temp;
    end
else
    for i = 1:length(ts)-1
        % Runge-Kutta 4
        % Looks like it is about 3-4x slower than euler forward. Should
        % compare the score when i get it working later.
        [f1,U(:,i+1),inputs(:,i+1)] = URED_APC(ts(i),         Y(:,i),                 ref,param);
        [f2,~,~]                    = URED_APC(ts(i)+1/3*dt,  Y(:,i)+dt/3*f1,         ref,param);
        [f3,~,~]                    = URED_APC(ts(i)+2/3*dt,  Y(:,i)+dt*(-f1/3+f2),   ref,param);
        [f4,~,~]                    = URED_APC(ts(i)+dt,      Y(:,i)+dt*(f1-f2+f3),   ref,param);
        Y(:,i+1)                    = Y(:,i) + dt/8*(f1+3*f2+3*f3+f4);
    end
end
toc

%% PLOTTING
figure
hold on
plot(Y(1,:),Y(2,:))
plot(ref{1}(ts),ref{2}(ts),':')
grid on
ylabel('X_1(t)')
xlabel('X_2(t)')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(1:2,:)-[ref{1}(ts);ref{2}(ts)]) ))
hold off

figure
hold on
plot(ts,Y(1,:)-ref{1}(ts))
plot(ts,Y(7,:),':')
grid on
xlabel('t')
legend('e1', 'ehat1')
title('Error estimation 1')
hold off

figure
hold on
plot(ts,Y(2,:)-ref{2}(ts))
plot(ts,Y(8,:),':')
grid on
xlabel('t')
legend('e2', 'ehat2')
title('Error estimation 2')
hold off

figure
hold on
plot(ts,Y(1,:))
plot(ts,ref{1}(ts),':')
grid on
xlabel('t')
legend('Actual path', 'Reference path')
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts),':')
grid on
xlabel('t')
legend('Actual path', 'Reference path')
hold off

figure
hold on
plot(ts,inputs(1,:))
grid on
xlabel('t')
title('F')
hold off

figure
hold on
plot(ts,inputs(2,:))
grid on
xlabel('t')
title('phidot')
hold off

figure
hold on
plot(ts,U(1,:))
plot(ts,U(2,:))
plot(ts,U(3,:),':')
plot(ts,U(4,:),':')
yline(param.ubar)
yline(-param.ubar)
grid on
xlabel('t')
title('U')
legend('u1','u2','u_{d,1}','u_{d,2}','u_{sat}')
hold off

%%
function [res,U,inputs] = URED_APC(t,y,ref,params)
    Y = y(1:6);
    ehat = y(7:12);
    rho = y(13:14);
   
    Y_ref = cellfun(@(c) c(t), ref).';
    % e = Y(1:2)-Y_ref(1:2);
    e = Y-Y_ref;

    % URED
    % sigma(1) = e
    % sigma(2) = edot
    % z(1)(t) = ehat
    % z(2)(t) = ehatdot
    % Thus sigma(2) is the estimator error.
    % Estimate velocity error
    sigmax = ehat(1)-e(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + params.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*params.mu*sigmax + 1.5*params.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(1:2:4)-params.k.*phix).';
    
    sigmay = ehat(2)-e(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + params.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*params.mu*sigmay + 1.5*params.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(2:2:4)-params.k.*phiy).';
    ehatdot1 = [ehatdotx;ehatdoty];
    ehatdot1 = ehatdot1(:);

    % extended observer.
    sigmax = ehat(3)-ehatdot1(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + params.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*params.mu*sigmax + 1.5*params.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(3:2:6)-params.k.*phix).';
    
    sigmay = ehat(4)-ehatdot1(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + params.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*params.mu*sigmay + 1.5*params.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(4:2:6)-params.k.*phiy).';
    ehatdot2 = [ehatdotx;ehatdoty];
    ehatdot2 = ehatdot2(:);
    
    ehatdot = [ehatdot1(1:2); ehatdot2];

    % Feedback control
    % 1. Linear filter
    s(1) = params.a*sat([e(1); ehat([3 5])], params.ebar);
    s(2) = params.a*sat([e(2); ehat([4 6])], params.ebar);

    % 2. APC

    xi = s.'./rho;
    ud = -params.K*Jf(xi).*Tf(xi);
    u = sat(ud, params.ubar);
    U = [u; ud];

    rhodot = (u-ud)./xi - params.lu.*(rho - params.rho_inf);

    ydot = params.A*Y + params.B*u;

    res = [ydot; ehatdot; rhodot];

    % Vechicle control signals.
    u1 = sqrt(Y(3)^2+Y(4)^2);
    phi = atan((Y(3)*Y(6)-Y(4)*Y(5))/(u1^3)*params.l);
    u1dot = (Y(3)*Y(5)+Y(4)*Y(6))/u1;

    phidot = (params.l*sqrt(Y(3)^2 + Y(4)^2) * (u(2) * Y(3) - Y(4) * u(1))...
        * (Y(3)^2 + Y(4)^2) - 3*(Y(5) * Y(3) + Y(6) * Y(4)) * (Y(6) * Y(3)...
        - Y(5) * Y(4))) / ((Y(3)^2 + Y(4)^2)^3 + params.l^2 * (Y(6) * Y(3)...
        - Y(5)*Y(4))^2);

    F = (params.M1+params.M2+params.M2*tan(phi)^2)*u1dot + ...
        params.M2*phidot*sin(phi)*u1/(cos(phi)^3);

    inputs = [F; phidot];
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

function res = RMS(X)
    % X is array of scalars or vectors
    temp = zeros(size(X,1),1);
    for i=1:size(X,1)
        temp(i) = X(:,i).'*X(:,i);
    end
    res = sqrt(1/length(X) * sum(temp) );
end