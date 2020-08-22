figure('Name','KPBC Analysis','NumberTitle','off');

axKPBC1 = subplot(2,3,1);
plot(time, hip_pos(I));
xlabel("Time (s)")
ylabel("Angle (deg)")
title("Hip Position")

axKPBC2 = subplot(2,3,2);
plot(time,knee_joint_vel(I).*U_pbc_k(I), time, ankle_joint_vel(I).*U_pbc_a(I));
xlabel("Time (s)")
ylabel("Power (J/s)")
legend("PBC Knee", "PBC Ankle");

axKPBC3 = subplot(2,3,3);
plot(time,U_pbc_k(I), time, U_pbc_a(I));
xlabel("Time (s)")
ylabel("Torque (N m)")
legend("PBC Knee", "PBC Ankle");

% ax118 = subplot(2,3,4);
% plot(time, hip_vel(I));
% xlabel("Time (s)")
% ylabel("Angle (deg/s)")
% title("Hip Velocity")
% linkaxes([ax66, ax118],'x')

axKPBC5 = subplot(2,3,5);
plot(time,Esys_int_out(I),time, Eref(I))
xlabel("Time (s)")
ylabel("Energy (Joules)")
legend("E_{sys} out","E_{ref}")

axKPBC6 = subplot(2,3,6);
plot(time, knee_joint_vel(I),time,ankle_joint_vel(I));
legend("Knee Vel","Ankle Vel")
xlabel("Time (s)")
ylabel("Angular Velocity (deg/s)")

linkaxes([axKPBC1,axKPBC2,axKPBC3,axKPBC5,axKPBC6],'x')

% ax171 = subplot(2,3,6);
% plot(time,PushOff(I))
% xlabel("Time (s)")
% title("PushOff")