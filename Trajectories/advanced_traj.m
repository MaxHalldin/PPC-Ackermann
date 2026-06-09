clear all
clc
% Data file for advanced trajectory.

tspan = [0; 46]; % Seconds
sampling_freq = 2e4; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

A = 2;

vias = [ 
  % x     y   t
    0     0   0;
    0     0   1;
    0     0   2;
    4*A   0   10;
    % 4.5*A 0   11;
    5*A   0   12;
    % 5*A .5*A  13;
    5*A  1*A  14;
    5*A  4*A  20;
    % 5*A 4.5*A 21;
    5*A  5*A  22;
    % 4.5*A 5*A 23;
    4*A  5*A  24;
    A    5*A  30;
    % .5*A 5*A  31;
    0    5*A  32;
    % 0   4.5*A 33;
    0    4*A  34;
    0     0   42;
    0     0   44;
    0     0   46;
];

traj = zeros(2,length(ts));
traj(1,:) = spline(vias(:,3),vias(:,1),ts);
traj(2,:) = spline(vias(:,3),vias(:,2),ts);

save(".\Trajectories\adv_traj.mat")

%%
figure
plot(traj(1,:),traj(2,:))
figure
plot(ts, traj(1,:))
figure
plot(ts,traj(2,:))

% y0 = [0; 0; 0; 0; 0; 0]