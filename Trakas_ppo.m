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
tspan = [0; 2.5*pi]; % Seconds

sampling_freq = 1e3; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

params = struct;
params.e_inf = 0.05;
params.lambda = 2;
params.r0 = 10;
params.r_inf = params.lambda(1)*params.e_inf;
params.gamma = [1.1; 1.0];
params.lo = 100;
params.A = [0 1; 0 0];

ehat_0 = [-0.2; 0.3];

%%

Y = zeros(length(ehat_0),length(ts));
Y(:,1) = ehat_0;
U = zeros(1,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

euler = true;
tic
tot = length(ts);
if euler
    for i = 1:length(ts)-1
        % This is the forward Euler method. 
        Y(:,i+1) = Y(:,i) + dt*edotppo(ts(i),Y(:,i),ref,params);
    end
else
    for i = 1:length(ts)-1
        % Runge-Kutta 4
        % Looks like it is about 3-4x slower than euler forward. Should
        % compare the score when i get it working later.
        f1      = edotppo(ts(i),         Y(:,i),                 ref,params);
        f2      = edotppo(ts(i)+1/3*dt,  Y(:,i)+dt/3*f1,         ref,params);
        f3      = edotppo(ts(i)+2/3*dt,  Y(:,i)+dt*(-f1/3+f2),   ref,params);
        f4      = edotppo(ts(i)+dt,      Y(:,i)+dt*(f1-f2+f3),   ref,params);
        Y(:,i+1)    = Y(:,i) + dt/8*(f1+3*f2+3*f3+f4);
    end
end
toc

%% PLOTTING
figure
hold on
plot(ts,Y(1,:))
plot(ts,ref{1}(ts))
grid on
xlabel('t')
legend('Estimated', 'Reference')
title(sprintf('y\nMSE = %0.3e', mean( (Y(1,:)-ref{1}(ts)).^2 )))
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts))
grid on
xlabel('t')
legend('Estimated', 'Reference')
title(sprintf('dy/dt\nMSE = %0.3e', mean( (Y(2,:)-ref{2}(ts)).^2 )))
hold off

%%
function res = edotppo(t,ehat,ref,pars)

    Y_ref = cellfun(@(c) c(t), ref).';
    r = (pars.r0 - pars.r_inf)*exp(-pars.lo*t) + pars.r_inf;

    % PPO
    e = Y_ref; %Assuming perfect controller.
    eta = ( e(1) - ehat(1) )/r;
    Gamma = pars.gamma./(r.^(1:2).');
    ehatdot = pars.A*ehat + Gamma*Tf(eta);

    res = ehatdot;
end

%% Helper functions

function res = Tf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    % assert(isbetween(x,-1,1), "Argument is out of bounds.");
    I = ones(size(x));
    res = 0.5 * log( I+x ./ (I-x) );
end