function [u_k, u_a, t_out, dt_out] = PDControlTest(...
    knee_pos, ankle_pos, knee_des, ankle_des, dt, kp_k, kd_k, kp_a, kd_a, time_in, filter_coeff)


%% Variables Defined
%Persistent variables used to store data between iterations
persistent knee_pos_prev ankle_pos_prev dknee_prev dankle_prev t t_prev

    if isempty(t_prev)
        t_prev = time_in;
    end
    if isempty(knee_pos_prev) || ((time_in - t_prev)>1)
        knee_pos_prev = 0;
        ankle_pos_prev = 0;
        t = 0;
        dknee_prev = 0;
        dankle_prev = 0;
    end


%% Initialization
    
    % Software/Hardware limits
%     ANKLE_POS_MIN_LIM  = -35;
%     ANKLE_POS_MAX_LIM = 35;
%     ANKLE_VEL_MIN_LIM = -200; %deg/s
%     ANKLE_VEL_MAX_LIM = 200; %deg/s
%     KNEE_POS_MIN_LIM = 2;
%     KNEE_POS_MAX_LIM = 85;
%     KNEE_VEL_MIN_LIM = -400; %deg/s
%     KNEE_VEL_MAX_LIM = 400; %deg/s

%% Control
    %Calculate joint velocities using weighted backwards difference
    t = dt + t;
    dknee = (1-filter_coeff)*dknee_prev + filter_coeff*(knee_pos-knee_pos_prev)/dt;
    dankle = (1-filter_coeff)*dankle_prev + filter_coeff*(ankle_pos-ankle_pos_prev)/dt;
       
    %[u_k,u_a] = PDControl(kp_k,kd_k,knee_des,knee_pos,dknee,kp_a,kd_a,ankle_des,ankle_pos,dankle);
    u_k = kp_k*(knee_des-knee_pos)   + kd_k*(-dknee);
    u_a = kp_a*(ankle_des-ankle_pos) + kd_a*(-dankle);
    
    max_torque = 60;
    u_k = min(max_torque, max(-max_torque, u_k));
    u_a = min(max_torque, max(-max_torque, u_a));
    
%% Storing persistent variables for next iteration
    knee_pos_prev = knee_pos;
    ankle_pos_prev = ankle_pos;
    t_out = t;
    t_prev = time_in;
    dknee_prev = dknee;
    dankle_prev = dankle;
    dt_out = dt;
end

% %% Helper functions
% function [Knee_torque_command,Ankle_torque_command] = PDControl(kp_k, kd_k, q_kstar, q_k, qdot_k, kp_a, kd_a, q_astar, q_a, qdot_a)
%     %This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero.
%     Knee_torque_command = kp_k*(q_kstar-q_k) + kd_k*(-qdot_k);
%     Ankle_torque_command = kp_a*(q_astar-q_a) + kd_a*(-qdot_a);
% end
% 
% function y = Saturate(x,x1,x2)
%     %Function to prevent the desired joint angles from changing to fast. 
%     %Works via saturation
%     y=min(x,max(x1,x2));
%     y=max(y,min(x1,x2));
% end
