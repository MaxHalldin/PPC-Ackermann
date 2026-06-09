clear all; 
clc; 
close all;

%% 

data = load("circ_traj.mat"); 
% data = load("eight_traj.mat"); 
% data = load("adv_traj.mat"); 
traj = data.traj;
ts = data.ts;
dt = data.dt;

param = struct;

param.l = 0.5;
param.M1 = 2;
param.M2 = 2;
param.Fbar = 10; % Placeholder
param.phidotbar = 10; % Placeholder

param.e_inf = 0.05*[1 1].';
param.lu = 1.5;
param.lambda = [2 1].';
param.a = flip(poly(-1*param.lambda)); %Row vector % Need to have the negative here for the roots(lambdas).
param.rho_0 = 5*[1 1].';
param.rho_inf = param.lambda(1)*param.e_inf;
param.K = 1.2;
param.ubar = 5;
param.ebar = 5;

param.A = [zeros(4,2) eye(4,4); zeros(2,6)];
param.B = [zeros(4,2); eye(2,2)];

param.mu = 1;
param.L = 2.5;
param.k = [2*sqrt(3);6]; % values borrowed from article.
cond1 = ( 0<param.k(1) && param.k(1)<2*sqrt(param.L) ) && ...
        ( param.k(2) > param.k(1)^2/4 + 4*param.L/(param.k(2)^2) );
cond2 = ( param.k(1)>2*sqrt(param.L) && param.k(2) > 2*param.L );
assert(cond1||cond2);

y0 = [-0.01; 1.99; 0.1; 0.01; 0; 0];

ehat_0 = [0 0 0 0 0 0].';
xer0 = [y0; ehat_0; param.rho_0];

%%

Y = zeros(length(xer0),length(ts));
Y(:,1) = xer0;
U = zeros(4,length(ts));
inputs = zeros(4,length(ts));

% assert(r_0 > abs(y_0(1)-yr(0)-ehat_0(1)))

tic
for i = 1:length(ts)-1
    % This is the forward Euler method. 
    [temp,U(:,i+1),inputs(:,i+1)] = URED_APC(ts(i), Y(:,i), traj(1:2,i) ,param);
    Y(:,i+1) = Y(:,i) + dt*temp;
end
toc

%% PLOTTING

figure
hold on
plot(Y(1,:),Y(2,:))
plot(traj(1,:),traj(2,:),':')
grid on
ylabel('X_1(t)')
xlabel('X_2(t)')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(1:2,:)-traj(1:2,:)) ))
hold off

figure
hold on
plot(ts,Y(1,:)-traj(1,:))
plot(ts,Y(7,:),':')
grid on
xlabel('t')
legend('e1', 'ehat1')
title('Error estimation 1')
hold off

figure
hold on
plot(ts,Y(2,:)-traj(2,:))
plot(ts,Y(8,:),':')
grid on
xlabel('t')
legend('e2', 'ehat2')
title('Error estimation 2')
hold off

figure
hold on
plot(ts,Y(1,:))
plot(ts,traj(1,:),':')
grid on
xlabel('t')
legend('Actual path', 'Reference path')
title("x-coordinate tracking")
hold off

figure
hold on
plot(ts,Y(2,:))
plot(ts,traj(2,:),':')
grid on
xlabel('t')
legend('Actual path', 'Reference path')
title("y-coordinate tracking")
hold off

figure
hold on
plot(ts,inputs(2,:))
plot(ts,inputs(1,:))
ylim([-1.2*param.Fbar 1.2*param.Fbar])
grid on
xlabel('t')
title('F')
hold off

figure
hold on
plot(ts,inputs(4,:))
plot(ts,inputs(3,:))
ylim([-1.2*param.phidotbar 1.2*param.phidotbar])
grid on
xlabel('t')
title('phidot')
hold off

figure
hold on
plot(ts,U(1,:))
plot(ts,U(2,:))
plot(ts,U(3,:),':')
plot(ts,U(4,:),':')
grid on
xlabel('t')
title('U')
legend('u1','u2','u_{d,1}','u_{d,2}')
hold off

figure
hold on
plot(ts,Y(13,:))
plot(ts,-Y(13,:))
grid on
xlabel('t')
title('rho_1')
hold off

figure
hold on
plot(ts,Y(14,:))
plot(ts,-Y(14,:))
grid on
xlabel('t')
title('rho_2')
hold off

%%
function [res,U,inputs] = URED_APC(~,y,ref,params)
    Y = y(1:6);
    ehat = y(7:12);
    rho = y(13:14);
   
    e = Y(1:2)-ref;

    % URED
    sigmax = ehat(1)-e(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + params.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*params.mu*sigmax + 1.5*params.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(1:2:4)-params.k.*phix).';
    
    sigmay = ehat(2)-e(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + params.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*params.mu*sigmay + 1.5*params.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(2:2:4)-params.k.*phiy).';
    ehatdot1 = [ehatdotx;ehatdoty];
    ehatdot1 = ehatdot1(:);

    % extended observer.
    sigmax = ehat(3)-ehatdot1(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + params.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*params.mu*sigmax + 1.5*params.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(3:2:6)-params.k.*phix).';
    
    sigmay = ehat(4)-ehatdot1(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + params.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*params.mu*sigmay + 1.5*params.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(4:2:6)-params.k.*phiy).';
    ehatdot2 = [ehatdotx;ehatdoty];
    ehatdot2 = ehatdot2(:);
    
    ehatdot = [ehatdot1(1:2); ehatdot2];

    % Feedback control
    % 1. Linear filter
    % s(1) = params.a*sat([e(1); ehat([3 5])], params.ebar);
    % s(2) = params.a*sat([e(2); ehat([4 6])], params.ebar);

    s(1) = params.a*[e(1); ehat([3 5])];
    s(2) = params.a*[e(2); ehat([4 6])];

    % 2. APC
    xi = s.'./rho;
    ud = -params.K*Jf(xi).*Tf(xi);

    % Vechicle control signals.
    
    u1  = sqrt(Y(3)^2+Y(4)^2);
    wd  = (Y(3)*Y(5)+Y(4)*Y(6))/u1;
    phi = atan(params.l * (Y(3)*Y(6)-Y(4)*Y(5)) / u1^3 );

    phidotd = (params.l*sqrt(Y(3)^2 + Y(4)^2) * (ud(2)*Y(3) - Y(4)*ud(1))...
        * (Y(3)^2 + Y(4)^2) - 3*(Y(5)*Y(3) + Y(6)*Y(4)) * (Y(6)*Y(3)...
        - Y(5)*Y(4))) / ((Y(3)^2 + Y(4)^2)^3 + params.l^2 * (Y(6) * Y(3)...
        - Y(5)*Y(4))^2);
    
    wdotd = 1/(sqrt(Y(3)^2+Y(4)^2))^3 * ...
        ( (Y(3)^2+Y(4)^2) * (Y(5)^2+Y(6)^2 + ud(1)*Y(4)+ud(2)*Y(3)) - ...
        ( Y(3)*Y(5) + Y(4)*Y(6) )^2 );


    Fd = (params.M1+params.M2+params.M2*tan(phi)^2)*wd + ...
        params.M2*phidotd*sin(phi)*u1/(cos(phi)^3);
    
    % Imput Saturation
    phidot = sat(phidotd, params.phidotbar);
    F = sat(Fd, params.Fbar);
    w = (F-params.M2*phidot*sin(phi)/(cos(phi)^3)) / ...
        (params.M1 + params.M2 + params.M1*tan(phi));

    a = (wdotd-w^2*u1^3)/(u1+eps) - w^2 - u1^4*tan(phi)/params.l^2;
    b = 3*u1^3*w*tan(phi) + u1^3*phidotd*(1 + tan(phi)^2)/params.l;

    u(1) = (a-b)/(2*Y(4)+eps);
    u(2) = (a+b)/(2*Y(3)+eps);
    u = u.';
    
    cond = ( (abs(u(1))<abs(ud(1))) && (abs(u(2))<abs(ud(2))) );
    u = cond*u + ~cond*ud; % Use ud if imputs are not saturated.
    U = [u; ud];

    inputs = [F; Fd; phidot; phidotd];

    rhodot = (u-ud)./xi - params.lu.*(rho - params.rho_inf);

    ydot = params.A*Y + params.B*u;

    res = [ydot; ehatdot; rhodot];

end

%% Helper functions

function res = Tf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    % assert(isbetween(x,-1,1), "Argument is out of bounds.");
    I = ones(size(x));
    res = 0.5 * log( I+x ./ (I-x) );
end

function res = Jf(x) % x in (-1,1)
    % x can be scalar or vector (or a matrix)
    % assert(isbetween(x,-1,1), "Argument is out of bounds.");
    I = ones(size(x));
    res = 1./(I-x.^2);
end

function res = sat(x, xbar)
    % x can be vector or scalar. If x is a scalar, xbar must also be. If x 
    % is a vector, then xbar can either be a scalar or a vector of same 
    % size as x.
    res = sign(x).*min(abs(x), xbar);
end

function res = RMS(X)
    % X is array of scalars or vectors
    temp = zeros(size(X,1),1);
    for i=1:size(X,1)
        temp(i) = X(:,i).'*X(:,i);
    end
    res = sqrt(1/length(X) * sum(temp) );
end