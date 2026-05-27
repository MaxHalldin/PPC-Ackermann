clear all; 
clc; 
close all;

%% 

data = load("circ_traj.mat"); 
ref = data.ref;

tspan = [0; 10]; % Seconds
sampling_freq = 1e4; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

param = struct;

param.dt = dt;

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

param.k = [6; 4.5; 4.182]; % values borrowed from article.
param.DGF = @(x)sign(x)*(abs(x)^(1/2)+abs(x)^(3/2));
param.DGFprime = @(x)(3*abs(x)+1)/(2*sqrt(abs(x)) + eps );

y0 = [-0.01; 1.99; 0.99; -0.01; 0; 0];

ehat_0 = [0 0 0 0 0 0].';
xer0 = [y0; ehat_0; param.rho_0];

%%

Y = zeros(length(xer0),length(ts));
Y(:,1) = xer0;
U = zeros(4,length(ts));
inputs = zeros(2,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

tic
for i = 1:length(ts)-1
    % This is the forward Euler method. 
    [temp,U(:,i+1),inputs(:,i+1)] = RED_APC(ts(i), Y(:,i), ref, param);
    Y(:,i+1) = Y(:,i) + dt*temp;
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
function [res,U,inputs] = RED_APC(t,y,ref,params)
    Y = y(1:6);
    ehat = y(7:12);
    rho = y(13:14);
   
    Y_ref = cellfun(@(c) c(t), ref);
    % e = Y(1:2)-Y_ref(1:2);
    e = Y-Y_ref;
    de = e(1:2)-ehat(1:2);

    % RED
    vx = [1/params.k(3)*params.DGF(params.k(3)^2*de(1));
         2*params.DGF(params.k(3)^2*de(1))*params.DGFprime(params.k(3)^2*de(1))];
    ehatdotx = ([ehat(3); 0] + params.k(1:2).*vx).';
    
    vy = [1/params.k(3)*params.DGF(params.k(3)^2*de(2));
         2*params.DGF(params.k(3)^2*de(2))*params.DGFprime(params.k(3)^2*de(2))];
    ehatdoty = ([ehat(4); 0] + params.k(1:2).*vy).';
    
    ehatdot1 = [ehatdotx;ehatdoty];
    ehatdot1 = ehatdot1(:);

    % Extending observer
    
    dedot = params.dt*ehatdot1(3:4);

    % This will need to be tuned experimentally
    k = [10; 10; 0.2].*params.k; 

    vx = [1/k(3)*params.DGF(k(3)^2*dedot(1));
         2*params.DGF(k(3)^2*dedot(1))*params.DGFprime(k(3)^2*dedot(1))];
    ehatdotx = ([ehat(5); 0] + k(1:2).*vx).';
    
    vy = [1/k(3)*params.DGF(k(3)^2*dedot(2));
         2*params.DGF(k(3)^2*dedot(2))*params.DGFprime(k(3)^2*dedot(2))];
    ehatdoty = ([ehat(6); 0] + k(1:2).*vy).';
    
    ehatdot2 = [ehatdotx;ehatdoty];
    ehatdot2 = ehatdot2(:);

    % Adding it all together
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
    % Need to implement the capping of u based on F and phidot


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

function res = sat_pos(x, xbar)
    % Assuming xbar is positive
    temp = min(abs(x), xbar);
    res = max(zeros(size(temp)), temp);
    
end

function res = RMS(X)
    % X is array of scalars or vectors
    temp = zeros(size(X,1),1);
    for i=1:size(X,1)
        temp(i) = X(:,i).'*X(:,i);
    end
    res = sqrt(1/length(X) * sum(temp) );
end