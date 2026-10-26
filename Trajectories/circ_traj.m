% Data file for circular trajectory

tspan = [0; 20]; % Seconds
sampling_freq = 1e4; %Hz
dt = 1/sampling_freq;
ts = tspan(1):dt:tspan(2);

A = 1;
c = 2*pi/tspan(2);

ref = cell(6,1);
ref{1} = @(t)A*sin(c*t);
ref{2} = @(t)A*cos(c*t);
ref{3} = @(t)A*c*cos(c*t);
ref{4} = @(t)A*-c*sin(c*t);
ref{5} = @(t)A*-c^2*sin(c*t);
ref{6} = @(t)A*-c^2*cos(c*t);

traj = cell2mat(cellfun(@(x) x(ts), ref, 'uni', 0 ));

save(".\Trajectories\circ_traj.mat")

%%
% plot(traj(1,:),traj(2,:))

% RED
% y0 = [-0.01; 0.99; 0.3; -0.01; 0; 0];
