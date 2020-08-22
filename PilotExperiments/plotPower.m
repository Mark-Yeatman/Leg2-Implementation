figure('Name','Power','NumberTitle','off');

ax11 = subplot(2,2,1);
plot(time,knee_joint_vel.*U_lin_spring_k, time, ankle_joint_vel.*U_lin_spring_a);
hold on;
plot(time,knee_joint_vel.*U_lin_damp_k, time, ankle_joint_vel.*U_lin_damp_a);
xlabel('Time (s)')
ylabel('Power (J/s)')
legend('Spring Knee', 'Spring Ankle', 'Damper Knee', 'Damper Ankle');

ax12 = subplot(2,2,2);
plot(time,knee_joint_vel.*U_pbc_k, time, ankle_joint_vel.*U_pbc_a);
xlabel('Time (s)')
ylabel('Power (J/s)')
legend('PBC Knee', 'PBC Ankle');

ax13 = subplot(2,2,3);
plot(time,knee_joint_vel.*U_stop_k, time, ankle_joint_vel.*U_stop_a);
xlabel('Time (s)')
ylabel('Power (J/s)')
legend('Knee Stop', 'Ankle Stop');

ax14 = subplot(2,2,4);
plot(time,knee_joint_vel.*Knee_torque_command, time, ankle_joint_vel.*Ankle_torque_command);
xlabel('Time (s)')
ylabel('Power (J/s)')
legend('Knee Net', 'Ankle Net');
linkaxes([ax11 ax12 ax13 ax14],'x');
axis tight
