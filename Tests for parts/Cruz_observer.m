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

sampling_freq = 1e5; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

params = struct;
params.A = [0 1; 0 0];
params.mu = 1;
params.L = 2.5;
params.k = [2*sqrt(3);6]; % values borrowed from article.
cond1 = ( 0<params.k(1) && params.k(1)<2*sqrt(params.L) ) && ...
        ( params.k(2) > params.k(1)^2/4 + 4*params.L/(params.k(2)^2) );
cond2 = ( params.k(1)>2*sqrt(params.L) && params.k(2) > 2*params.L );
assert(cond1||cond2);

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
        Y(:,i+1) = Y(:,i) + dt*edotURED(ts(i),Y(:,i),ref,params);
        % clc; fprintf("%0.1e",i/(length(ts)-1))
    end
else
    for i = 1:length(ts)-1
        % Runge-Kutta 4
        % Looks like it is about 3-4x slower than euler forward. Should
        % compare the score when i get it working later.
        f1      = edotURED(ts(i),         Y(:,i),                 ref,params);
        f2      = edotURED(ts(i)+1/3*dt,  Y(:,i)+dt/3*f1,         ref,params);
        f3      = edotURED(ts(i)+2/3*dt,  Y(:,i)+dt*(-f1/3+f2),   ref,params);
        f4      = edotURED(ts(i)+dt,      Y(:,i)+dt*(f1-f2+f3),   ref,params);
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
legend('Estimated','Reference')
title(sprintf('dt/dt\nMSE = %0.3e', mean( (Y(1,:)-ref{1}(ts)).^2 )))
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts))
grid on
xlabel('t')
legend('Estimated','Reference')
title(sprintf('dt/dt\nMSE = %0.3e', mean( (Y(2,:)-ref{2}(ts)).^2 )))
hold off

%%
function res = edotURED(t,ehat,ref,param)
    % URED
    % sigma(1) = y(t)
    % sigma(2) = ydot(t)
    % z(1)(t) is the esimator of y(t)
    % z(2)(t) is the esimator of ydot(t)
    % Thus sigma(2) is the estimator error.
    
    e = ref{1}(t); % I might be studpid but 
    sigma = ehat(1)-e(1);
    
    phi = [ abs(sigma)^(1/2)*sign(sigma) + param.mu*abs(sigma)^(3/2)*sign(sigma);
            0.5*sign(sigma) + 2*param.mu*sigma + 1.5*param.mu^2*abs(sigma)^2*sign(sigma) ];
    
    ehatdot = param.A*ehat-param.k.*phi;
    res = ehatdot;
end