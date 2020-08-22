figure('Name','Load Cell','NumberTitle','off');

% ax90 = subplot(3,1,1);
% plot(time, load_cell_x_force(I),time, load_cell_y_force(I), time, load_cell_z_force(I),...
%      time, load_cell_x_moment(I),  time, load_cell_y_moment(I),  time, load_cell_z_moment(I));
% legend("F_x","F_y","F_z", "M_x", "M_y", "M_z")
% xlabel("Time (s)")
% ylabel("Force and Moment (N and N m)")

ax91 = subplot(3,1,1);
Fnorm = vecnorm([load_cell_x_force,load_cell_y_force,load_cell_z_force]')';
plot(time, load_cell_x_force(I),time, load_cell_y_force(I), time, load_cell_z_force(I), time, Fnorm(I))
legend("F_x","F_y","F_z", "F")
xlabel("Time (s)")
ylabel("Force (N)")

ax92 = subplot(3,1,2);
Mnorm = vecnorm([load_cell_x_moment,load_cell_y_moment,load_cell_z_moment]')';
plot(time, load_cell_x_moment(I),  time, load_cell_y_moment(I),  time, load_cell_z_moment(I), time, Mnorm(I));
legend("M_x", "M_y", "M_z", "M")
xlabel("Time (s)")
ylabel("Moment (N m)")
linkaxes([ax91,ax92],'x')

ax171 = subplot(3,1,3);
plot(time,foot_contact(I))
xlabel("Time (s)")
title("Foot Contact")