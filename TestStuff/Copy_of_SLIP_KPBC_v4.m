function [Knee_torque_command, Ankle_torque_command, deltaL, hip_pos, t_out, dt]...
    = SLIP_KPBC_v1(IMU_pitch, Knee_motor_position, Knee_joint_position, Ankle_motor_position, ...
    Ankle_joint_position, Iteration, Iteration_time, Knee_torque_sensor, Ankle_torque_sensor, ...
    Load_cell_x_force, Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,...
    Load_cell_z_moment, kp_knee, kd_knee, kp_ankle,kd_ankle, ankle_ind_of_hip,...
    knee_ind_of_hip, ankle_des_in, knee_des_in, time_in, filter_coeff, IMU_filter_coeff,...
    q_h_0, q_h_min, c, s_po, ...
    FC, lf, la, ls, lt, k, d, L0, hip_vel, knee_vel, ankle_vel, M, Eref)

%
% Inputs in addition to sensors:

% kp_ankle,kd_ankle,kp_knee,kd_knee are joint values

% ankle_ind_of_hip and knee_ind_of_hip are boolean values (true or false),
% --can be checkboxes in LabView with default TRUE.

% ankle_des_cte and knee_des_cte are desired angles for ankle and knee,
% controlled from LabView.


%% Variables Defined
%Persistent variables used to store data between iterations
persistent hip_pos_prev knee_pos_prev ankle_pos_prev dknee_prev dankle_prev dhip_prev ...
           t t_prev IMU_pitch_prev 

    if isempty(t_prev)
        t_prev = time_in;
    end
    if isempty(knee_pos_prev) || ((time_in - t_prev)>1)
        hip_pos_prev = 0;
        knee_pos_prev = 0;
        ankle_pos_prev = 0;
        t = 0;
        dknee_prev = 0;
        dankle_prev = 0;
        dhip_prev = 0;
        IMU_pitch_prev = 0;
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
    knee_pos = - Knee_motor_position; %reverse knee sign convention of biomechanics versus biped modeling
    ankle_pos = Ankle_motor_position;

    x = zeros(6,1);
    x(1) =  hip_pos;
    x(2) =  knee_pos;         
    x(3) =  ankle_pos;
    x(4) =  hip_vel;
    x(5) =  -knee_vel; %because its not calculated using a finite difference method, its just given to us right now, change for implementation
    x(6) =  ankle_vel;
    
    %Calculate joint velocities usindg weighted backwards difference
    dt = Iteration_time;
    t = Iteration_time + t;
    dhip = (1-filter_coeff)*dhip_prev + filter_coeff*(hip_pos-hip_pos_prev)/dt;
    dknee = (1-filter_coeff)*dknee_prev + filter_coeff*(knee_pos-knee_pos_prev)/dt;
    dankle = (1-filter_coeff)*dankle_prev + filter_coeff*(ankle_pos-ankle_pos_prev)/dt;
    
    if FC %Use SLIP Embedding Controller
        %[Knee_torque_command,Ankle_torque_command ,deltaL] = SLIPControl(ankle_pos, knee_pos, hip_pos, lf, la, ls, lt, k, L0);
        %Virtual linear spring from hip to foot
        fc = lf;
        %x should be 
        %l1 is shank length
        %l2 is thigh length
        %lf = 0; %at heel


        J = [0;(-1).*lt.*(2.*fc.*cosd(x(2)+x(3))+2.*ls.*sind(x(2))+la.*sind(x(2)+x(3))).*(4.*fc.^2+la.^2+4.*ls.^2+4.*lt.^2+8.*ls.*lt.*cosd(x(2))+4.*la.*ls.*cosd(x(3))+4.*la.*lt.*cosd(x(2)+x(3))+(-8).*fc.*ls.*sind(x(3))+(-8).*fc.*lt.*sind(x(2)+x(3))).^(-1/2);(-1).*(4.*fc.^2+la.^2+4.*ls.^2+4.*lt.^2+8.*ls.*lt.*cosd(x(2))+4.*la.*ls.*cosd(x(3))+4.*la.*lt.*cosd(x(2)+x(3))+(-8).*fc.*ls.*sind(x(3))+(-8).*fc.*lt.*sind(x(2)+x(3))).^(-1/2).*(2.*fc.*ls.*cosd(x(3))+2.*fc.*lt.*cosd(x(2)+x(3))+la.*ls.*sind(x(3))+la.*lt.*sind(x(2)+x(3)))];
        L  = (fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cosd(x(2))+la.*ls.*cosd(x(3))+la.*lt.*cosd(x(2)+x(3))+(-2).*fc.*ls.*sind(x(3))+(-2).*fc.*lt.*sind(x(2)+x(3))).^(1/2);
        Ldot = (1/2).*(fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(2))+la.*ls.*cos(x(3))+la.*lt.*cos(x(2)+x(3))+(-2).*fc.*ls.*sin(x(3))+(-2).*fc.*lt.*sin(x(2)+x(3))).^(-1/2).*((-2).*ls.*lt.*sin(x(2)).*x(5)+(-2).*fc.*ls.*cos(x(3)).*x(6)+(-1).*la.*ls.*sin(x(3)).*x(6)+(-2).*fc.*lt.*cos(x(2)+x(3)).*(x(5)+x(6))+(-1).*la.*lt.*sin(x(2)+x(3)).*(x(5)+x(6)));
        deltaL = L - L0;
        u = -k.*(deltaL).*J - d.*(Ldot).*J;
        Espring = 1/2*k*deltaL^2;
        Knee_torque_command  = u(2);
        Ankle_torque_command = u(3);

        if true %energy shaping
            Ehip = HipEnergy(x(1:3),x(4:6),M);
            u_KPBC = KPBC(x(4:6),[1;1], Ehip, Eref);
        end
        
    else %Use PD Controller       
        %Calculate torque output using PD controller
        [Knee_torque_command,Ankle_torque_command] = PDControl(...
            kp_knee,kd_knee,knee_des_in,knee_pos,dknee, ...
            kp_ankle,kd_ankle,ankle_des_in,ankle_pos,dankle);
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

function [Knee_torque_command,Ankle_torque_command,deltaL,E] = SLIPControl(ankle, knee, hip, fc, la, ls, lt, ks, L0)
    %Virtual linear spring from hip to foot
    %x should be 
    %l1 is shank length
    %l2 is thigh length
    %lf = 0; %at heel
    x = zeros(3,1);
    x(1) = ankle+knee-hip;
    x(2) = knee;
    x(3) = ankle;

    J = [0;(-1).*lt.*(2.*fc.*cosd(x(2)+x(3))+2.*ls.*sind(x(2))+la.*sind(x(2)+x(3))).*(4.*fc.^2+la.^2+4.*ls.^2+4.*lt.^2+8.*ls.*lt.*cosd(x(2))+4.*la.*ls.*cosd(x(3))+4.*la.*lt.*cosd(x(2)+x(3))+(-8).*fc.*ls.*sind(x(3))+(-8).*fc.*lt.*sind(x(2)+x(3))).^(-1/2);(-1).*(4.*fc.^2+la.^2+4.*ls.^2+4.*lt.^2+8.*ls.*lt.*cosd(x(2))+4.*la.*ls.*cosd(x(3))+4.*la.*lt.*cosd(x(2)+x(3))+(-8).*fc.*ls.*sind(x(3))+(-8).*fc.*lt.*sind(x(2)+x(3))).^(-1/2).*(2.*fc.*ls.*cosd(x(3))+2.*fc.*lt.*cosd(x(2)+x(3))+la.*ls.*sind(x(3))+la.*lt.*sind(x(2)+x(3)))];
    L  = (fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cosd(x(2))+la.*ls.*cosd(x(3))+la.*lt.*cosd(x(2)+x(3))+(-2).*fc.*ls.*sind(x(3))+(-2).*fc.*lt.*sind(x(2)+x(3))).^(1/2);
    deltaL = L - L0;
    u = -ks*(deltaL)*J;
    E = 1/2*ks*deltaL^2;
    
    Knee_torque_command  = u(2);
    Ankle_torque_command = u(3);
end

function E = ProsEnergy(q_a, q_k, qdot_a, qdot_k, lf, la, ls, Mf, Ma, Ms)
    %Only defined for stance phase where the prothesis is constrained to the ground
%     x(1) = ankle+knee-hip;
%     x(2) = knee;
%     x(3) = ankle;
%     
%     KE = [(1/2).*(x(6).*((la.*Ms.*(la+ls.*cosd(x(3))+lt.*cosd(x(2)+x(3)))+(1/4).*lf.*Ms.*(lf+(-2).*ls.*sind(x(3))+(-2).*lt.*sind(x(2)+x(3)))).*x(4)+(la.*Ms.*(la+ls.*cosd(x(3)))+(1/4).*lf.*Ms.*(lf+(-2).*ls.*sind(x(3)))).*x(5)+(la.^2.*Ms+(1/4).*lf.^2.*Ms).*x(6))+x(5).*(((1/2).*ls.*Ms.*((1/2).*ls+lt.*cosd(x(2)))+Ms.*(la+ls.*cosd(x(3))).*(la+ls.*cosd(x(3))+lt.*cosd(x(2)+x(3)))+(1/4).*Ms.*(lf+(-2).*ls.*sind(x(3))).*(lf+(-2).*ls.*sind(x(3))+(-2).*lt.*sind(x(2)+x(3)))).*x(4)+((1/4).*ls.^2.*Ms+Ms.*(la+ls.*cosd(x(3))).^2+(1/4).*Ms.*(lf+(-2).*ls.*sind(x(3))).^2).*x(5)+(la.*Ms.*(la+ls.*cosd(x(3)))+(1/4).*lf.*Ms.*(lf+(-2).*ls.*sind(x(3)))).*x(6))+x(4).*(((1/4).*lt.^2.*Mt+Ms.*((1/2).*ls+lt.*cosd(x(2))).^2+Ms.*(la+ls.*cosd(x(3))+lt.*cosd(x(2)+x(3))).^2+lt.^2.*Ms.*sind(x(2)).^2+(1/4).*Ms.*(lf+(-2).*ls.*sind(x(3))+(-2).*lt.*sind(x(2)+x(3))).^2).*x(4)+((1/2).*ls.*Ms.*((1/2).*ls+lt.*cosd(x(2)))+Ms.*(la+ls.*cosd(x(3))).*(la+ls.*cosd(x(3))+lt.*cosd(x(2)+x(3)))+(1/4).*Ms.*(lf+(-2).*ls.*sind(x(3))).*(lf+(-2).*ls.*sind(x(3))+(-2).*lt.*sind(x(2)+x(3)))).*x(5)+(la.*Ms.*(la+ls.*cosd(x(3))+lt.*cosd(x(2)+x(3)))+(1/4).*lf.*Ms.*(lf+(-2).*ls.*sind(x(3))+(-2).*lt.*sind(x(2)+x(3)))).*x(6)))];
%     PE = (-1/2).*g.*(lt.*(2.*Ms+Mt).*cosd(x(1))+ls.*Ms.*cosd(x(1)+x(2)));
end

function E = HipEnergy(x,v,M)
   ToeP = [lf.*cosd(x(1)+x(2)+x(3))+lt.*sind(x(1))+ls.*sind(x(1)+x(2))+la.*sind(x(1)+x(2)+x(3)),(-1).*lt.*cosd(x(1))+(-1).*ls.*cosd(x(1)+x(2))+(-1).*la.*cosd(x(1)+x(2)+x(3))+lf.*sind(x(1)+x(2)+x(3)),0,0,0,sind(x(1)+x(2)+x(3))]; 
   HipP = -ToeP; %because of the switch from HipFrame to ToeFrame
   ToeV = [-(v(1)+v(2)+v(3))*lf.*sind(x(1)+x(2)+x(3)) + v(1)*lt.*cosd(x(1)) + (v(1)+v(2))*ls.*cosd(x(1)+x(2)) + (v(1)+v(2)+v(3))*la.*cosd(x(1)+x(2)+x(3)),...
       v(1).*lt.*sind(x(1)) + (v(1)+v(2)).*ls.*sind(x(1)+x(2)) + (v(1)+v(2)+v(3)).*la.*sind(x(1)+x(2)+x(3)) + (v(1)+v(2)+v(3))*lf.*cosd(x(1)+x(2)+x(3)),...
       0,...
       0,...
       0,...
       (v(1)+v(2)+v(3))*cosd(x(1)+x(2)+x(3))]; 
   HipV = -ToeV; %because of the switch from HipFrame to ToeFrame
   
   E = norm(HipV(1:2),2)*M*1/2 + HipP(2)*M*9.81;
end

function u_KPBC = KPBC(v,omega,E,Eref)
    %Kinetic-Passivity Based Control
    u_KPBC = -omega*v(E-Eref);
end

function y = Saturate(x,x1,x2)
    %Function to prevent the desired joint angles from changing to fast. 
    %Works via saturation
    y=min(x,max(x1,x2));
    y=max(y,min(x1,x2));
end
