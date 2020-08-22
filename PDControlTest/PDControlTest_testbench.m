%This is a test bench for the matlab function. This script is
%necessary in the creation of the C code. 

%Initialize all variables used in the control
knee_pos = single(0.06);
ankle_pos = single(-10);
knee_des = single(0);
ankle_des = single(0);
dt = single(0.0001);
kp_k = single(1);
kd_k = single(0);
kp_a = single(12);
kd_a = single(3);
time_in = single(1);
filter_coeff = single(0.9);

%Call the control function.
[u_k, u_a, t_out, dt_out] = PDControlTest(...
    knee_pos, ankle_pos, knee_des, ankle_des, dt, kp_k, kd_k, kp_a, kd_a, time_in, filter_coeff)