clear all
clc
% Data file for advanced trajectory.

tspan = [0; 92]; % Seconds
sampling_freq = 5e3; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

A = 2;

vias = [ 
  % x     y   t
    0     0   0;
    0     0   2;
    0     0   4;
    4*A   0   20;
    5*A   0   24;
    5*A  1*A  28;
    5*A  4*A  40;
    5*A  5*A  44;
    4*A  5*A  48;
    A    5*A  60;
    0    5*A  64;
    0    4*A  68;
    0     0   84;
    0     0   88;
    0     0   92;
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