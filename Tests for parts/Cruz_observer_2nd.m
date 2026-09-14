clear all; 
clc; 
close all;

%% I suspect that my problem is that I need to sample faster than ode45 
% does for integrating
params = struct;

A = 1;
c = 1;
y1t  = @(t)A*sin(c*t);
y2t  = @(t)A*cos(c*t);
y1tp = @(t)A*c*cos(c*t);
y2tp = @(t)A*-c*sin(c*t);
y1tb = @(t)A*-c^2*sin(c*t);
y2tb = @(t)A*-c^2*cos(c*t);
ref = {y1t y2t y1tp y2tp y1tb y2tb};
tspan = [0; 7]; % Seconds

sampling_freq = 1e4; %Hz
params.dt = 1/sampling_freq;
ts = tspan(1):params.dt:tspan(2);

params.A = [0 1; 0 0];
params.mu = 1;
params.L = 2.5;
params.k = [2*sqrt(3);6]; % values borrowed from article.
cond1 = ( 0<params.k(1) && params.k(1)<2*sqrt(params.L) ) && ...
        ( params.k(2) > params.k(1)^2/4 + 4*params.L/(params.k(2)^2) );
cond2 = ( params.k(1)>2*sqrt(params.L) && params.k(2) > 2*params.L );
assert(cond1||cond2);

ehat_0 = [-0.2; 0.3; 0.1; 0.1; 0.1; 0.1];

%%

Y = zeros(length(ehat_0),length(ts));
% dYdt = zeros(size(Y));
Y(:,1) = ehat_0;
U = zeros(1,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

tic
for i = 1:length(ts)-1
    Y(:,i+1) = Y(:,i) + params.dt*edotURED(ts(i),Y(:,i),{ref{1:2}},params);
end
toc

%% PLOTTING
figure
hold on
plot(ts,Y(1,:))
plot(ts,ref{1}(ts),':')
grid on
ylabel('X_1(t)')
xlabel('X_2(t)')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(1,:)-ref{1}(ts)) ))
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts),':')
grid on
ylabel('X(t)')
xlabel('t')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(2,:)-ref{2}(ts)) ))
hold off

figure
hold on
plot(Y(1,:),Y(2,:))
plot(ref{1}(ts),ref{2}(ts),':')
grid on
ylabel('Y(t)')
xlabel('t')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS( Y(1:2,:)-[ref{1}(ts);ref{2}(ts)] ) ))
hold off

figure
hold on
plot(ts,Y(1,:))
plot(ts,ref{1}(ts),':')
grid on
ylabel('X_1(t)')
xlabel('X_2(t)')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(1,:)-ref{1}(ts)) ))
hold off

figure
hold on
plot(ts,Y(3,:))
plot(ts,ref{3}(ts),':')
grid on
ylabel('t')
xlabel('dx')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(3,:)-ref{3}(ts)) ))
hold off

figure
hold on
plot(ts,Y(4,:))
plot(ts,ref{4}(ts),':')
grid on
ylabel('t')
xlabel('dy')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(4,:)-ref{4}(ts)) ))
hold off

figure
hold on
plot(ts,Y(5,:))
plot(ts,ref{5}(ts),':')
grid on
ylabel('t')
xlabel('dx')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(5,:)-ref{5}(ts)) ))
hold off

figure
hold on
plot(ts,Y(6,:))
plot(ts,ref{6}(ts),':')
grid on
ylabel('t')
xlabel('dy')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(6,:)-ref{6}(ts)) ))
hold off

%%
function ehatdot = edotURED(t,y,ref,param)
    
    ehat = y(1:6);
   
    Y_ref = cellfun(@(c) c(t), ref).';
    e = Y_ref;

    % URED
    sigmax = ehat(1)-e(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + param.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*param.mu*sigmax + 1.5*param.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(1:2:4)-param.k.*phix).';
    
    sigmay = ehat(2)-e(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + param.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*param.mu*sigmay + 1.5*param.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(2:2:4)-param.k.*phiy).';
    ehatdot1 = [ehatdotx;ehatdoty];
    ehatdot1 = ehatdot1(:);

    % extended observer.
    sigmax = ehat(3)-ehatdot1(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + param.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*param.mu*sigmax + 1.5*param.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(3:2:6)-param.k.*phix).';
    
    sigmay = ehat(4)-ehatdot1(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + param.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*param.mu*sigmay + 1.5*param.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(4:2:6)-param.k.*phiy).';
    ehatdot2 = [ehatdotx;ehatdoty];
    ehatdot2 = ehatdot2(:);
    
    ehatdot = [ehatdot1(1:2); ehatdot2];
end

function res = RMS(X)
    % X is array of scalars or vectors
    temp = zeros(size(X,1),1);
    for i=1:size(X,1)
        temp(i) = X(:,i).'*X(:,i);
    end
    res = sqrt(1/length(X) * sum(temp) );
end