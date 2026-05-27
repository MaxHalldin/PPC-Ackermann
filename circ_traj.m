% Data file for circular trajectory

A = 2;
c = pi/5;

ref = cell(6,1);
ref{1}  = @(t)A*sin(c*t);
ref{2}  = @(t)A*cos(c*t);
ref{3} = @(t)A*c*cos(c*t);
ref{4} = @(t)A*-c*sin(c*t);
ref{5} = @(t)A*-c^2*sin(c*t);
ref{6} = @(t)A*-c^2*cos(c*t);

save("circ_traj.mat")