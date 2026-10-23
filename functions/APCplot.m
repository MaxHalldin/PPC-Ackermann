function APCplot(traj, ts, Y, U, inputs, S,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    traj
    ts
    Y
    U
    inputs
    S
    param
end

arguments (Output)
end

figure
hold on
plot(Y(1,:),Y(2,:))
plot(traj(1,:),traj(2,:),':')
grid on
ylabel('X_1(t)')
xlabel('X_2(t)')
legend('Estimated path', 'True path')
title(sprintf('Path tracking\nRMS = %0.3e', RMS(Y(1:2,:)-[traj(1,:);traj(2,:)]) ))
hold off

figure
subplot(2,1,1);
plot(ts,Y(1,:)-traj(1,:),ts,Y(7,:),':')
lgd1 = legend('$e_1$','$\hat{e}_1$');
xlabel('t')
grid on
title(sprintf('RMSE_1 = %0.3e', RMS( Y(7,:)-(Y(1,:)-traj(1,:))) ))
subplot(2,1,2);
plot(ts,Y(2,:)-traj(2,:),ts,Y(8,:),':')
lgd2 = legend('$e_2$','$\hat{e}_2$');
xlabel('t')
grid on
title(sprintf('RMSE_2 = %0.3e', RMS( Y(8,:)-(Y(2,:)-traj(2,:))) ))
sgtitle('Error estimation')
set(lgd1, 'interpreter', 'latex')
set(lgd2, 'interpreter', 'latex')

figure
subplot(2,1,1);
plot(ts,Y(1,:),ts,traj(1,:),':')
grid on
xlabel('t')
ylabel('y_1')
legend('Actual path', 'Reference path')
title(sprintf('RMS_1 = %0.3e', RMS(Y(1,:)-traj(1,:))) )
subplot(2,1,2);
plot(ts,Y(2,:),ts,traj(2,:),':')
grid on
xlabel('t')
ylabel('y_2')
legend('Actual path', 'Reference path')
title(sprintf('RMS_2 = %0.3e', RMS(Y(2,:)-traj(2,:))) )
sgtitle("Coordinate tracking")

figure
subplot(2,1,1);
hold on
plot(ts,inputs(2,:))
plot(ts,inputs(1,:))
yline([-param.Fbar,param.Fbar],'--')
ylim([-3*param.Fbar 3*param.Fbar])
grid on
xlabel('t')
title('F')
legend('Unconstrained', 'Constrained')
hold off
subplot(2,1,2);
hold on
plot(ts,inputs(4,:))
plot(ts,inputs(3,:))
yline([-param.phidotbar,param.phidotbar],'--')
ylim([-3*param.phidotbar 3*param.phidotbar])
grid on
xlabel('t')
t = '$\dot{\phi}$';
title(t,'interpreter', 'latex')
legend('Unconstrained', 'Constrained')
hold off
sgtitle('Input saturation')

figure
subplot(2,2,1);
hold on
plot(ts,S(1,:))
plot(ts,Y(1,:)-traj(1,:))
plot(ts,Y(13,:),'b')
plot(ts,-Y(13,:),'b')
grid on
xlabel('t')
legend('s_1','e_1','\rho_1')
hold off
subplot(2,2,3);
hold on
plot(ts,S(2,:))
plot(ts,Y(2,:)-traj(2,:))
plot(ts,Y(14,:),'b')
plot(ts,-Y(14,:),'b')
grid on
xlabel('t')
legend('s_2','e_2','\rho_2')
hold off
subplot(2,2,2);
start = 20000;
hold on
plot(ts(start:end),S(1,start:end))
plot(ts(start:end),Y(1,start:end)-traj(1,start:end))
plot(ts(start:end),Y(13,start:end),'b')
plot(ts(start:end),-Y(13,start:end),'b')
grid on
xlabel('t')
legend('s_1','e_1','\rho_1')
hold off
subplot(2,2,4);
hold on
plot(ts(start:end),S(2,start:end))
plot(ts(start:end),Y(2,start:end)-traj(2,start:end))
plot(ts(start:end),Y(14,start:end),'b')
plot(ts(start:end),-Y(14,start:end),'b')
grid on
xlabel('t')
legend('s_2','e_2','\rho_2')
hold off
sgtitle('Adaptive performance funnels')

figure
subplot(2,1,1);
hold on
plot(ts,U(1,:))
plot(ts,U(3,:))
grid on
xlabel('t')
legend('u_1','u_{d,1}')
hold off
subplot(2,1,2);
hold on
plot(ts,U(2,:))
plot(ts,U(4,:))
grid on
xlabel('t')
legend('u_2','u_{d,2}')
hold off
sgtitle('Flat control inputs')

end