load robustFCtestData.mat

[t,I] = sort(robustFCtest.ankle.t_out);

ankle_pos = robustFCtest.ankle.JointPosition(I);
ankle_vel = robustFCtest.ankle.JointVelocity(I);
knee_pos = robustFCtest.knee.JointPosition(I);
knee_vel = robustFCtest.knee.JointVelocity(I);
t = t - min(t);

%magic cutoff from inspection
I = t > 56.25 & t < 110;
t = t(I);
ankle_pos = ankle_pos(I);
knee_pos = knee_pos(I);
ankle_vel = ankle_vel(I);
knee_vel = knee_vel(I);

figure
f1 = subplot(2,2,1);
plot(t,ankle_pos);
xlabel('t (seconds)')
ylabel('ankle joint angle')
title('ankle position')

f2 = subplot(2,2,2);
plot(t,knee_pos);
xlabel('t (seconds)')
ylabel('knee joint angle')
title('knee position')

f3 = subplot(2,2,3);
plot(t,ankle_vel)
xlabel('t (seconds)')
ylabel('ankle angular vel')
title('ankle velocity')

f4 = subplot(2,2,4);
plot(t,knee_vel)
xlabel('t (seconds)')
ylabel('knee angular vel')
title('knee velocity')

IMU_pitch=single(0);
Knee_joint_position=single(4.39); 
Ankle_joint_position=single(1); 
time_in=single(1);
Iteration_time=single(1);
Load_cell_x_force=single(1);
Load_cell_y_force=single(1); 
Load_cell_z_force=single(1); 
Load_cell_x_moment=single(1); 
Load_cell_y_moment=single(1);
Load_cell_z_moment=single(1);
kp_knee=single(-1);
kd_knee=single(0);
kp_ankle=single(-1);
kd_ankle=single(0);
ankle_des_in=single(0);
knee_des_in=single(0);
vel_filter_coeff=single(1);
IMU_filter_coeff=single(1);
SLIP_ON = single(0);
lf=single(0);
la=single(1); 
ls=single(1); 
lt=single(1); 
k=single(2); 
d=single(0); 
L0=single(1);
KPBC_ON = single(0);
pbc_gain_ankle = single(0); 
pbc_gain_knee = single(0);
M = single(0); 
Eref = single(0);
knee_stop_low = single(2);
knee_stop_high = single(105);
ankle_stop_low = single(-35);
ankle_stop_high = single(35);
%matlab GUI tool butterworth
for i = 1:length(ankle_pos)
[Knee_torque_command(i), Ankle_torque_command(i), deltaL, hip_pos, t_out(i), dt, Esys, Esys_int_out]...
= SLIP_KPBC(IMU_pitch, knee_pos(i), ankle_pos(i), ...
            t(i), Iteration_time, ...
            Load_cell_x_force, Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment,...
            kp_knee, kd_knee, kp_ankle, kd_ankle, ankle_des_in, knee_des_in,...
            vel_filter_coeff, IMU_filter_coeff,...
            SLIP_ON, lf, la, ls, lt, k, d, L0, ...
            KPBC_ON, pbc_gain_ankle, pbc_gain_knee, M, Eref,...
            knee_stop_low, knee_stop_high, ankle_stop_low, ankle_stop_high );
end
ankle_pos_filt_mat = Ankle_torque_command;
knee_pos_filt_mat = Knee_torque_command;
ankle_vel_filt_mat = gradient(ankle_pos_filt_mat,t);
knee_vel_filt_mat = gradient(knee_pos_filt_mat,t);
subplot(f1);
hold on
plot(t,ankle_pos_filt_mat)
legend('raw','mat')

subplot(f2);
hold on
plot(t,knee_pos_filt_mat)
legend('raw','mat')

%custom exponential weighting
vel_filter_coeff = 0.5;
dt = t(2)-t(1);
ankle_vel_filt_exp = zeros(size(ankle_pos));
knee_vel_filt_exp = zeros(size(ankle_pos));
for i = 2:length(ankle_pos)
    ankle_vel_filt_exp(i) = (1-vel_filter_coeff)*ankle_vel_filt_exp(i-1) + vel_filter_coeff*(ankle_pos(i)-ankle_pos(i-1))/dt;
    knee_vel_filt_exp(i) = (1-vel_filter_coeff)*knee_vel_filt_exp(i-1) + vel_filter_coeff*(knee_pos(i)-knee_pos(i-1))/dt;
end

%custom butterworth 2nd order
NZEROS = 2;
NPOLES = 2;
GAIN = 7.485478157e+01;
xv = zeros(NZEROS+1,2);
yv = zeros(NPOLES+1,2);

dt = t(2)-t(1);
ankle_vel_filt_butt = zeros(size(ankle_pos));
knee_vel_filt_butt = zeros(size(ankle_pos));
for i = 3:length(ankle_pos)
    xv(1,:) = xv(2,:); 
    xv(2,:) = xv(3,:); 
    xv(3,:) = [(ankle_pos(i)-ankle_pos(i-1))/dt;(knee_pos(i)-knee_pos(i-1))/dt]./ GAIN;
    yv(1,:) = yv(2,:); 
    yv(2,:) = yv(3,:); 
    yv(3,:) =  (xv(1,:) + xv(3,:)) + 2 * xv(2,:)+ ( -0.7008967812 * yv(1,:)) + (  1.6474599811 * yv(2,:));
    ankle_vel_filt_butt(i) = yv(3,1);
    knee_vel_filt_butt(i) = yv(3,2);
end

subplot(f3);
hold on
plot(t,ankle_vel_filt_exp,t,ankle_vel_filt_butt,t,ankle_vel_filt_mat)
legend('raw','exp','butt','mat')

subplot(f4);
hold on
plot(t,knee_vel_filt_exp,t,knee_vel_filt_butt,t,knee_vel_filt_mat)
legend('raw','exp','butt','mat')
