%File constructing the trajectory for the ackermann robot.
clearvars -except traj; clc; close all

function res = ref(t,c,A)
    res = cell(6,1);
    if t < pi/(2*c)
        res{1} = A*(1-cos(c*t));
        res{2} = A*(0+sin(c*t));
        res{3} = A*c*sin(c*t);
        res{4} = A*c*cos(c*t);
        res{5} = A*c^2*cos(c*t);
        res{6} = -A*c^2*sin(c*t);
        return
    elseif t < pi/(c) 
        res{1} = A*(1-cos(c*t));
        res{2} = A*(2-sin(c*t));
        res{3} = A*c*sin(c*t);
        res{4} = -A*c*cos(c*t);
        res{5} = A*c^2*cos(c*t);
        res{6} = -A*c^2*sin(c*t);
        return
    elseif t < 2*pi/(c) 
        res{1} = A*(3+cos(c*t));
        res{2} = A*(2-sin(c*t));
        res{3} = A*c*sin(c*t);
        res{4} = -A*c*cos(c*t);
        res{5} = -A*c^2*cos(c*t);
        res{6} = A*c*2*sin(c*t);
        return
    elseif t < 4*pi/(c)
        res{1} = 4*A;
        res{2} = A*(4-c*t/pi);
        res{3} = 0;
        res{4} = -A*c/pi;
        res{5} = 0;
        res{6} = 0;
        return
    else
        res{1} = A*(3+cos(c*t));
        res{2} = A*(0-sin(c*t));
        res{3} = -A*c*sin(c*t);
        res{4} = -A*c*cos(c*t);
        res{5} = -A*c^2*cos(c*t);
        res{6} = -A*c^2*sin(c*t);
    end
end

function res = y1(t,c,A)
    res = cell(6,1);
    if t < pi/(2*c)
        res{1} = A*(1-cos(c*t));
        res{2} = A*(0+sin(c*t));
        res{3} = A*c*sin(c*t);
        res{4} = A*c*cos(c*t);
        res{5} = A*c^2*cos(c*t);
        res{6} = -A*c^2*sin(c*t);
        return
    elseif t < pi/(c) 
        res{1} = A*(1-cos(c*t));
        res{2} = A*(2-sin(c*t));
        res{3} = A*c*sin(c*t);
        res{4} = -A*c*cos(c*t);
        res{5} = A*c^2*cos(c*t);
        res{6} = -A*c^2*sin(c*t);
        return
    elseif t < 2*pi/(c) 
        res{1} = A*(3+cos(c*t));
        res{2} = A*(2-sin(c*t));
        res{3} = A*c*sin(c*t);
        res{4} = -A*c*cos(c*t);
        res{5} = -A*c^2*cos(c*t);
        res{6} = A*c*2*sin(c*t);
        return
    elseif t < 4*pi/(c)
        res{1} = 4*A;
        res{2} = A*(4-c*t/pi);
        res{3} = 0;
        res{4} = -A*c/pi;
        res{5} = 0;
        res{6} = 0;
        return
    else
        res{1} = A*(3+cos(c*t));
        res{2} = A*(0-sin(c*t));
        res{3} = -A*c*sin(c*t);
        res{4} = -A*c*cos(c*t);
        res{5} = -A*c^2*cos(c*t);
        res{6} = -A*c^2*sin(c*t);
    end
end

c = pi/10;
A = 1;
t = 0:0.1:16*pi;

Y = cell(6,length(t));
for i = 1:length(t)
    Y(:,i) = ref(t(i),c,A);
end
traj = cell2mat(Y);

R = arrayfun(@(z) ref(z,c,A),t, UniformOutput=false);

% plot(traj(1,:),traj(2,:))
plot(t,traj(4,:))