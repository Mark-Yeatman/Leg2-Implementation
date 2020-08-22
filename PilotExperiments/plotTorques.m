figure('Name','Torques','NumberTitle','off');

ax1 = subplot(2,2,1);
plot(time,U_lin_spring_k(I), time, U_lin_spring_a(I));
hold on;
plot(time,U_lin_damp_k(I), time, U_lin_damp_a(I));
xlabel("Time (s)")
ylabel("Torque (N m)")
legend("Spring Knee", "Spring Ankle", "Damper Knee", "Damper Ankle");

ax2 = subplot(2,2,2);
plot(time,U_pbc_k(I), time, U_pbc_a(I));
xlabel("Time (s)")
ylabel("Torque (N m)")
legend("PBC Knee", "PBC Ankle");

ax3 = subplot(2,2,3);
plot(time,U_stop_k(I), time, U_stop_a(I));
xlabel("Time (s)")
ylabel("Torque (N m)")
legend("Knee Stop", "Ankle Stop");

ax4 = subplot(2,2,4);
plot(time,Knee_torque_command(I), time, Ankle_torque_command(I));
xlabel("Time (s)")
ylabel("Torque (N m)")
legend("Knee Net", "Ankle Net");
linkaxes([ax1 ax2 ax3 ax4],'x');
axis tight