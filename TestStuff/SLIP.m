function [Knee_torque_command, Ankle_torque_command, deltaL, hip_pos, t_out, dt]...
    = SLIP(IMU_pitch, Knee_motor_position, Knee_joint_position, Ankle_motor_position, ...
    Ankle_joint_position, Iteration, Iteration_time, Knee_torque_sensor, Ankle_torque_sensor, ...
    Load_cell_x_force, Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,...
    Load_cell_z_moment, kp_knee, kd_knee, kp_ankle,kd_ankle, ankle_ind_of_hip,...
    knee_ind_of_hip, ankle_des_in, knee_des_in, time_in, filter_coeff, IMU_filter_coeff,...
    q_h_0, q_h_min, c, s_po, ...
    FC, lf, la, ls, lt, k, d, L0, SLIP_ON)

%
% Inputs in addition to sensors:

% kp_ankle,kd_ankle,kp_knee,kd_knee are joint values

% ankle_ind_of_hip and knee_ind_of_hip are boolean values (true or false),
% --can be checkboxes in LabView with default TRUE.

% ankle_des_cte and knee_des_cte are desired angles for ankle and knee,
% controlled from LabView.


%% Variables Defined
%Persistent variables used to store data between iterations
persistent knee_pos_prev ankle_pos_prev dknee_prev dankle_prev ...
           t t_prev IMU_pitch_prev 

    if isempty(t_prev)
        t_prev = time_in;
    end
    if isempty(knee_pos_prev) || ((time_in - t_prev)>1)
        knee_pos_prev = 0;
        ankle_pos_prev = 0;
        t = 0;
        dknee_prev = 0;
        dankle_prev = 0;
        IMU_pitch_prev = 5;
    end


    %%
    %Initialization
    %Knee_torque_command = 0; 
    %Ankle_torque_command = 0; 
    deltaL = 0; 
    %hip_pos = 0; 
    %t_out = 0; 
    %dt = 0;
    
    % Software/Hardware limits
%     ANKLE_POS_MIN_LIM  = -35;
%     ANKLE_POS_MAX_LIM = 35;
%     ANKLE_VEL_MIN_LIM = -200; %deg/s
%     ANKLE_VEL_MAX_LIM = 200; %deg/s
%     KNEE_POS_MIN_LIM = 2;
%     KNEE_POS_MAX_LIM = 85;
%     KNEE_VEL_MIN_LIM = -400; %deg/s
%     KNEE_VEL_MAX_LIM = 400; %deg/s

    %Assign joint positions
    hip_pos = (1-IMU_filter_coeff)*IMU_pitch_prev + IMU_filter_coeff*IMU_pitch;
    knee_pos = Knee_motor_position;
    ankle_pos = Ankle_motor_position;

    %Calculate joint velocities usindg weighted backwards difference
    dt = Iteration_time;
    t = Iteration_time + t;

    dknee = (1-filter_coeff)*dknee_prev + filter_coeff*(knee_pos-knee_pos_prev)/dt;
    dankle = (1-filter_coeff)*dankle_prev + filter_coeff*(ankle_pos-ankle_pos_prev)/dt;
    
    if SLIP_ON && all( [lf,la,lt,ls] > 0.001)%Use SLIP Embedding Controller
        %reverse ankle and shift because sign convention of biomechanics versus biped modeling
        x = zeros(10);
        x(3) = hip_pos - knee_pos + ankle_pos;
        x(4) = -ankle_pos;
        x(5) = knee_pos;
        
        x(9) = -dankle;
        x(10) = dknee;
        
        J = [0;0;0;(lf.*ls.*cosd(x(4))+lf.*lt.*cosd(x(4)+x(5))+(-1).*la.*(ls.*sind(x(4))+lt.*sind(x(4)+x(5)))).*(la.^2+lf.^2+ls.^2+lt.^2+2.*(la.*ls.*cosd(x(4))+ls.*lt.*cosd(x(5))+la.*lt.*cosd(x(4)+x(5))+lf.*ls.*sind(x(4))+lf.*lt.*sind(x(4)+x(5)))).^(-1/2);lt.*(lf.*cosd(x(4)+x(5))+(-1).*ls.*sind(x(5))+(-1).*la.*sind(x(4)+x(5))).*(la.^2+lf.^2+ls.^2+lt.^2+2.*(la.*ls.*cosd(x(4))+ls.*lt.*cosd(x(5))+la.*lt.*cosd(x(4)+x(5))+lf.*ls.*sind(x(4))+lf.*lt.*sind(x(4)+x(5)))).^(-1/2);0;0;0];
        L = (la.^2+lf.^2+ls.^2+lt.^2+2.*(la.*ls.*cosd(x(4))+ls.*lt.*cosd(x(5))+la.*lt.*cosd(x(4)+x(5))+lf.*ls.*sind(x(4))+lf.*lt.*sind(x(4)+x(5)))).^(1/2);
        Ldot=(la.^2+lf.^2+ls.^2+lt.^2+2.*(la.*ls.*cosd(x(4))+ls.*lt.*cosd(x(5))+la.*lt.*cosd(x(4)+x(5))+lf.*ls.*sind(x(4))+lf.*lt.*sind(x(4)+x(5)))).^(-1/2).*((lf.*ls.*cosd(x(4))+lf.*lt.*cosd(x(4)+x(5))+(-1).*la.*(ls.*sind(x(4))+lt.*sind(x(4)+x(5)))).*x(9)+lt.*(lf.*cosd(x(4)+x(5))+(-1).*ls.*sind(x(5))+(-1).*la.*sind(x(4)+x(5))).*x(10));
        
        Fs = -k*(L-L0);
        Fd = -d*Ldot;
        F = Fs + Fd;
        max_force = 100;
        F = min(max_force, max(-max_force, F));
        
        u = J*F;
        
        Ankle_torque_command = -u(4);
        Knee_torque_command = u(5);
        
               
    else %Use PD Controller       
        %Calculate torque output usindddg PD controller
        [Knee_torque_command,Ankle_torque_command] = PDControl(...
            kp_knee,kd_knee,knee_des_in,knee_pos,dknee, ...
            kp_ankle,kd_ankle,ankle_des_in,ankle_pos,dankle);
        max_torque = 150;
        Knee_torque_command = min(max_torque, max(-max_torque, Knee_torque_command));
        Ankle_torque_command = min(max_torque, max(-max_torque, Ankle_torque_command));
    end

    %%
    %Storing persistent variables for next iteration
    knee_pos_prev = knee_pos;
    ankle_pos_prev = ankle_pos;
    t_out=t;
    t_prev=time_in;
    dknee_prev=dknee;
    dankle_prev=dankle;
    IMU_pitch_prev=hip_pos;
end
%function

%% Helper functions
function [Knee_torque_command,Ankle_torque_command] = PDControl(kp_k, kd_k, q_kstar, q_k, qdot_k, kp_a, kd_a, q_astar, q_a, qdot_a)
    %This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero.
    Knee_torque_command = kp_k*(q_kstar-q_k) + kd_k*(-qdot_k);
    Ankle_torque_command = kp_a*(q_astar-q_a) + kd_a*(-qdot_a);
end



