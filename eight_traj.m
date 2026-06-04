% Data file for circular trajectory

A = 1;
c = 10;

ref = cell(6,1);
ref{1}  = @(t)A*sin(c*t/10);
ref{2}  = @(t)A*sin(c*t/20);
ref{3} = @(t)A*c/10*cos(c*t/10);
ref{4} = @(t)A*c/20*cos(c*t/20);
ref{5} = @(t)A*-c^2/100*sin(c*t/10);
ref{6} = @(t)A*-c^2/400*sin(c*t/20);

save("eight_traj.mat")