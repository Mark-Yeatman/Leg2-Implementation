%This is a test bench (tb) for the matlab function. This script is
%necessary in the creation of the C code. 

%Initialized all variables used in the function below
IMU_pitch=single(1);
Knee_motor_position=single(1);
Knee_joint_position=single(1); 
Ankle_motor_position=single(1);
Ankle_joint_position=single(1); 
Iteration=single(1);
Iteration_time=single(1);
Knee_torque_sensor=single(1); 
Ankle_torque_sensor=single(1); 
Load_cell_x_force=single(1);
Load_cell_y_force=single(1); 
Load_cell_z_force=single(1); 
Load_cell_x_moment=single(1); 
Load_cell_y_moment=single(1);
Load_cell_z_moment=single(1);
kp_knee=single(1);
kd_knee=single(1);
kp_ankle=single(1);
kd_ankle=single(1);
ankle_ind_of_hip=single(1);
knee_ind_of_hip=single(1);
ankle_des_cte=single(1);
knee_des_cte=single(1);
time=single(1);
filter_coeff=single(1);
IMU_filter_coeff=single(1);
q_h_0=single(1);
q_h_min=single(1);
c=single(1);
s_po=single(1);
FC=single(1);


%Calls the control function. The control function houses all of the control
%logic and mathematics.

[Knee_torque_command, Ankle_torque_command, knee_des, ankle_des,t_out, hip,...
    ankle_d_a,ankle_d_d,dt]...
    = Holonomic_Controller_v1(IMU_pitch, Knee_motor_position,...
    Knee_joint_position, Ankle_motor_position, Ankle_joint_position, Iteration,...
    Iteration_time, Knee_torque_sensor, Ankle_torque_sensor, Load_cell_x_force,...
    Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,...
    Load_cell_z_moment,kp_knee,kd_knee,kp_ankle,kd_ankle,ankle_ind_of_hip,...
    knee_ind_of_hip,ankle_des_cte,knee_des_cte,time,filter_coeff,IMU_filter_coeff,...
   q_h_0,q_h_min,c,s_po,FC);



