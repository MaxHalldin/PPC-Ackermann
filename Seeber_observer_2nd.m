clear all; 
clc; 
close all;

data = load("circ_traj.mat"); 
ref = data.ref;
tspan = [0; 20]; % Seconds

params = struct;

sampling_freq = 1e4; %Hz
params.dt = 1/sampling_freq;
ts = tspan(1):params.dt:tspan(2);

params.A = [0 1; 0 0];
params.k = [6; 4.5; 4.182]; % values borrowed from article.
params.DGF = @(x)sign(x)*(abs(x)^(1/2)+abs(x)^(3/2));
params.DGFprime = @(x)(3*abs(x)+1)/(2*sqrt(abs(x)));

ehat_0 = [-0.2; 0.3; 0.1; 0.1; 0.1; 0.1];

%% Sim

Y = zeros(length(ehat_0),length(ts));
Y(:,1) = ehat_0;
U = zeros(1,length(ts));

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
xlabel('t')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(1,:)-ref{1}(ts)) ))
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,ref{2}(ts),':')
grid on
ylabel('X_2(t)')
xlabel('t')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(2,:)-ref{2}(ts)) ))
hold off

figure
hold on
plot(Y(1,:),Y(2,:))
plot(ref{1}(ts),ref{2}(ts),':')
grid on
ylabel('X_2(t)')
xlabel('X_1(t)')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS( Y(1:2,:)-[ref{1}(ts);ref{2}(ts)] ) ))
hold off

figure
hold on
plot(ts,Y(3,:))
plot(ts,ref{3}(ts),':')
grid on
xlabel('t')
ylabel('dx')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(3,:)-ref{3}(ts)) ))
hold off

figure
hold on
plot(ts,Y(4,:))
plot(ts,ref{4}(ts),':')
grid on
xlabel('t')
ylabel('dy')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(4,:)-ref{4}(ts)) ))
hold off

figure
hold on
plot(Y(3,:),Y(4,:))
plot(ref{3}(ts),ref{4}(ts),':')
grid on
ylabel('dX_2(t)')
xlabel('dX_1(t)')
legend('Estimated path', 'True path')
title(sprintf('Velocity tracking\nRMS = %0.3e', RMS( Y(3:4,:)-[ref{3}(ts);ref{4}(ts)] ) ))
hold off

figure
hold on
plot(ts,Y(5,:))
plot(ts,ref{5}(ts),':')
grid on
xlabel('t')
ylabel('d^2x')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(5,:)-ref{5}(ts)) ))
hold off

figure
hold on
plot(ts,Y(6,:))
plot(ts,ref{6}(ts),':')
grid on
xlabel('t')
ylabel('d^2y')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(6,:)-ref{6}(ts)) ))
hold off

figure
hold on
plot(Y(5,:),Y(6,:))
plot(ref{5}(ts),ref{6}(ts),':')
grid on
ylabel('d^2X_2(t)')
xlabel('d^2X_1(t)')
legend('Estimated path', 'True path')
title(sprintf('Acceleration tracking\nRMS = %0.3e', RMS( Y(5:6,:)-[ref{5}(ts);ref{6}(ts)] ) ))
hold off

figure
hold on
plot(Y(5,sampling_freq*4:end),Y(6,sampling_freq*4:end))
plot(ref{5}(ts(sampling_freq*4:end)),ref{6}(ts(sampling_freq*4:end)),':')
grid on
ylabel('d^2X_2(t)')
xlabel('d^2X_1(t)')
title(sprintf('Acceleration tracking\nRMS = %0.3e', RMS( Y(5:6,sampling_freq*4:end)-[ref{5}(ts(sampling_freq*4:end));ref{6}(ts(sampling_freq*4:end))] ) ))
legend('Estimated path', 'True path')
hold off

%% main
function Yhatdot = edotURED(t,ehat,ref,param)

    Yhat = ehat;
    Y_ref = cellfun(@(c) c(t), ref).';
    
    % RED
    e = Y_ref(1:2)-Yhat(1:2);
    vx = [1/param.k(3)*param.DGF(param.k(3)^2*e(1));
         2*param.DGF(param.k(3)^2*e(1))*param.DGFprime(param.k(3)^2*e(1))];
    Yhatdotx = (param.A*Yhat([1 3]) + param.k(1:2).*vx).';
    
    vy = [1/param.k(3)*param.DGF(param.k(3)^2*e(2));
         2*param.DGF(param.k(3)^2*e(2))*param.DGFprime(param.k(3)^2*e(2))];
    Yhatdoty = (param.A*Yhat([2 4]) + param.k(1:2).*vy).';
    
    Yhatdot1 = [Yhatdotx;Yhatdoty];
    Yhatdot1 = Yhatdot1(:);

    % Extending observer
    edot = param.dt*Yhatdot1(3:4);
    % This will need to be tuned experimentally
    k = [10; 10; 0.2].*param.k; 

    vx = [1/k(3)*param.DGF(k(3)^2*edot(1));
         2*param.DGF(k(3)^2*edot(1))*param.DGFprime(k(3)^2*edot(1))];
    Yhatdotx = (param.A*Yhat([3 5]) + k(1:2).*vx).';
    
    vy = [1/k(3)*param.DGF(k(3)^2*edot(2));
         2*param.DGF(k(3)^2*edot(2))*param.DGFprime(k(3)^2*edot(2))];
    Yhatdoty = (param.A*Yhat([4 6]) + k(1:2).*vy).';
    
    Yhatdot2 = [Yhatdotx;Yhatdoty];
    Yhatdot2 = Yhatdot2(:);

    % Adding it all together
    Yhatdot = [Yhatdot1(1:2); Yhatdot2];
end

function res = RMS(X)
    % X is array of scalars or vectors
    temp = zeros(size(X,1),1);
    for i=1:size(X,1)
        temp(i) = X(:,i).'*X(:,i);
    end
    res = sqrt(1/length(X) * sum(temp) );
end