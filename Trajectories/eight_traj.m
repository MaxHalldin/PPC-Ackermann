% Data file for circular trajectory

tspan = [0; 40*pi]; % Seconds
sampling_freq = 1e4; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

A = 2;
c = 1;

ref = cell(6,1);
ref{1}  = @(t)A*sin(c*t/10);
ref{2}  = @(t)A*sin(c*t/20);
ref{3} = @(t)A*c/10*cos(c*t/10);
ref{4} = @(t)A*c/20*cos(c*t/20);
ref{5} = @(t)A*-c^2/100*sin(c*t/10);
ref{6} = @(t)A*-c^2/400*sin(c*t/20);

traj = cell2mat(cellfun(@(x) x(ts), ref, 'uni', 0 ));

save(".\Trajectories\eight_traj.mat")

% plot(traj(1,:),traj(2,:))

% RED
% y0 = [-.01; .09; .3; -.01; 0; 0];