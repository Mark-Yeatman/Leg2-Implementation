%This is a test bench (tb) for the matlab function. This script is
%necessary in the creation of the C code. 

%Initialized all variables used in the function below
IMU_pitch=single(0);
Knee_joint_position=single(9.28); 
Ankle_joint_position=single(-3.61); 
time_in=single(1);
Iteration_time=single(1);
Load_cell_x_force=single(1);
Load_cell_y_force=single(1); 
Load_cell_z_force=single(1); 
Load_cell_x_moment=single(1); 
Load_cell_y_moment=single(1);
Load_cell_z_moment=single(1);
kp_knee=single(12);
kd_knee=single(0.3);
kp_ankle=single(40);
kd_ankle=single(0.25);
ankle_des_in=single(1-5);
knee_des_in=single(3);
vel_filter_coeff=single(1);
KPBC_filter_coeff=single(1);
SLIP_ON = single(1);
lf=single(0);
la=single(1); 
ls=single(1); 
lt=single(1); 
k=single(100); 
d=single(0); 
L0=single(0.8984);
KPBC_ON = single(0);
pbc_gain_ankle = single(0); 
pbc_gain_knee = single(0);
M = single(0); 
Eref = single(0);
knee_stop_low = single(2);
knee_stop_high = single(105);
ankle_stop_low = single(-35);
ankle_stop_high = single(35);
max_torque = single(100);
k_tilde = single(0);
Ignore_PushOff = single(10);
hip_thresh_sw = single(1);
F_thresh =single(10);
knee_des_minus = single(-2);
knee_des_plus = single(0);
ankle_des_minus = single(-2);
ankle_des_plus = single(0);
Command_State = single(0);
KPBC_max_torque = single(0);
KPBC_joint_level = single(0);
%Calls the control function. The control function houses all of the control
%logic and mathematics.

[Knee_torque_command, Ankle_torque_command, deltaL, hip_pos, Esys, Esys_integrate_out,...
        U_LIN_SPRING_K, U_LIN_SPRING_A, U_LIN_DAMP_K, U_LIN_DAMP_A, U_STOP_K, U_STOP_A, U_PBC_K, U_PBC_A, ...
        knee_des_out,ankle_des_out, foot_contact, stance, phase_out,IMU_LIVE_OUT]...
= SLIP_KPBC(IMU_pitch, Knee_joint_position, Ankle_joint_position, ...
            time_in, Iteration_time, ...
            Load_cell_x_force, Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment,...
            kp_knee, kd_knee, kp_ankle, kd_ankle, ankle_des_in, knee_des_in,...
            vel_filter_coeff, KPBC_filter_coeff,...
            SLIP_ON, lt, k, d, L0, ...
            KPBC_ON, pbc_gain_ankle, pbc_gain_knee, M, Eref,...
            knee_stop_low, knee_stop_high, ankle_stop_low, ankle_stop_high,...
            max_torque,...
            k_tilde,Ignore_PushOff, F_thresh,Command_State, KPBC_max_torque, KPBC_joint_level)



