clear all; 
clc; 
close all;

addpath("functions","Trajectories\");

% traj_path = ".\Trajectories\circ_traj";
% traj_path = ".\Trajectories\eight_traj";
traj_path = ".\Trajectories\adv_traj";

run(strcat(traj_path,".m"))
data = load(strcat(traj_path,".mat")); 

traj    = data.traj;
ts      = data.ts;

param = struct;
param.dt = data.dt;

% Physical properties
param.l = 0.5;
param.M1 = 10;
param.M2 = 10;
param.Fbar = 1.0; % Placeholder
param.phidotbar = 1.0; % Placeholder

% PPC properties
param.e_inf = 0.1*[1 1].';
param.lu = 1.2;
param.lambda = [2 1].';
param.a = flip(poly(-1*param.lambda)); %Row vector % Need to have the negative here for the roots(lambdas).
param.rho_0 = 5*[1 1].';
param.rho_inf = param.lambda(1)*param.e_inf;
param.K = 1.0;
param.ebar = 1;

% Model specific
param.A = [zeros(4,2) eye(4,4); zeros(2,6)];
param.B = [zeros(4,2); eye(2,2)];

% URED specific
param.mu = 1;
param.L = 2.5;
param.k = [2*sqrt(3);6]; % values borrowed from article.
cond1 = ( 0<param.k(1) && param.k(1)<2*sqrt(param.L) ) && ...
        ( param.k(2) > param.k(1)^2/4 + 4*param.L/(param.k(2)^2) );
cond2 = ( param.k(1)>2*sqrt(param.L) && param.k(2) > 2*param.L );
assert(cond1||cond2);

% Initial conditions
y0 = [-0.01; 0.1; 0.1; 0.01; 0; 0];
ehat_0 = [0 0 0 0 0 0].'; % Assuming that our initial guess is prefect
xer0 = [y0; ehat_0; param.rho_0];

%% Simulation

tot = length(ts);
Y = zeros(length(xer0),tot);
Y(:,1) = xer0;
U = zeros(4,tot);
inputs = zeros(4,tot);
S = zeros(2,tot);

tic
for i = 1:length(ts)-1
    % Standard forward Euler. 
    [temp,U(:,i+1),inputs(:,i+1),S(:,i+1)] = URED_APC(ts(i), Y(:,i), traj(1:2,i) ,param);
    Y(:,i+1) = Y(:,i) + param.dt*temp;
    % disp([num2str(100*i/tot,'%.2f'), ' %'])
end
toc

%% PLOTTING
APCplot(traj,ts,Y,U,inputs,S,param)