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
params.A = [0 1; 0 0];
params.k = [6; 4.5; 4.182]; % values borrowed from article.
params.DGF = @(x)sign(x)*(abs(x)^(1/2)+abs(x)^(3/2));
params.DGFprime = @(x)(3*abs(x)+1)/(2*sqrt(abs(x)));


ehat_0 = [-0.2; 0.3];

%%

Y = zeros(length(ehat_0),length(ts));
Y(:,1) = ehat_0;
U = zeros(1,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

euler = false;
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
legend('Estimated','Reference (w/ noise)')
title(sprintf('f\nMSE = %0.3e', RMS( Y(1,:)-ref{1}(ts) )))
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts))
grid on
xlabel('t')
legend('Estimated','Reference (w/ noise)')
title(sprintf('dt/dt\nMSE = %0.3e', RMS( Y(2,:)-ref{2}(ts) )))
hold off

%% main
function res = edotURED(t,ehat,ref,param)
    % RED
    yhat = ehat;
    e = ref{1}(t)-yhat(1);
    
    v = [1/param.k(3)*param.DGF(param.k(3)^2*e);
         2*param.DGF(param.k(3)^2*e)*param.DGFprime(param.k(3)^2*e)];

    yhatdot = param.A*yhat + param.k(1:2).*v;

    res = yhatdot;
end

function res = RMS(X)
    % X is array of scalars or vectors
    temp = zeros(size(X,1),1);
    for i=1:size(X,1)
        temp(i) = X(:,i).'*X(:,i);
    end
    res = sqrt(1/length(X) * sum(temp) );
end