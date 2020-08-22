function [Knee_torque_command, Ankle_torque_command, knee_des, ankle_des, ...
    t_out, hip,...
    ankle_d_a,ankle_d_d,dt]...
    = Holonomic_Controller_v1(IMU_pitch, Knee_motor_position,...
    Knee_joint_position, Ankle_motor_position, Ankle_joint_position, Iteration,...
    Iteration_time, Knee_torque_sensor, Ankle_torque_sensor, Load_cell_x_force,...
    Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,...
    Load_cell_z_moment,kp_knee,kd_knee,kp_ankle,kd_ankle,ankle_ind_of_hip,...
    knee_ind_of_hip,ankle_des_cte,knee_des_cte,time,filter_coeff,IMU_filter_coeff,...
   q_h_0,q_h_min,c,s_po,FC)
%
% Inputs in addition to sensors:

% kp_ankle,kd_ankle,kp_knee,kd_knee are joint values

% ankle_ind_of_hip and knee_ind_of_hip are boolean values (true or false),
% --can be checkboxes in LabView with default TRUE.

% ankle_des_cte and knee_des_cte are desired angles for ankle and knee,
% controlled from LabView.


%% Variables Defined
%Persistent variables used to store data between iterations
persistent knee_pos_p ankle_pos_p ankle_des_p knee_des_p knee_p ankle_p t time_p dknee_p dankle_p IMU_pitch_p q_h_m s_m state

if isempty(time_p)
    time_p=time;
end
if isempty(knee_pos_p) || ((time-time_p)>1)
    knee_pos_p=0;
    ankle_pos_p=0;
    ankle_des_p=0;
    knee_des_p=10;
    knee_p=10;
    ankle_p=0;
    t=0;
    dknee_p=0;
    dankle_p=0;
    IMU_pitch_p=5;
    q_h_m=-10.99;
    s_m=c;
    state=1;
end


%% Winters Data Arrays
a_knee=[-21.9472500000000,4.24525000000000,21.6720000000000,0;20.9295000000000,-6.97650000000000,7.71900000000000,0;15.9362500000000,-0.998250000000000,18.6640000000000,25.8830000000000;-59.8412500000000,13.6422500000000,64.8630000000000,0;-19.4405000000000,-1.40250000000000,38.4100000000000,-47.2960000000000;45.1070000000000,-7.15300000000000,0.456000000000000,0;1.42525000000000,0.0657500000000000,2.21000000000000,3.24500000000000];
a_ankle=[5.50050000000000,-0.880500000000000,-4.60000000000000,0;-8.98500000000000,3.77020000000000,5.26600000000000,4.65120000000000;0.812500000000000,-0.475000000000000,7.74100000000000,2.81250000000000;-2.19975000000000,0.320750000000000,9.62000000000000,0;-14.7575000000000,0.878500000000001,-0.745000000000000,-24.2440000000000;23.0932000000000,-3.91420000000000,-19.9240000000000,0;-29.8020000000000,9.93400000000000,-0.0560000000000000,0;0.703500000000000,-0.234500000000000,-0.525000000000000,0;-2.60850000000000,0.869500000000000,1.21400000000000,0;-0.729000000000000,0.0190000000000000,0.580000000000000,-1.34400000000000];
knee_ind=[1 139 400 530 720 848 976 1001];
ankle_ind=[1 61 205 330 440 550 653 846 898 977 1001];

% q_h_0=19.33;
% q_h_min=-10.99;
% c=.53;

%%
% To be used as inputs later on

ankle_pos_min_lim=-35;
ankle_pos_max_lim=35;
ankle_vel_min_lim=-200; %deg/s
ankle_vel_max_lim=200; %deg/s
knee_pos_min_lim=2;
knee_pos_max_lim=85;
knee_vel_min_lim=-400; %deg/s
knee_vel_max_lim=400; %deg/s


%% 

hip=(1-IMU_filter_coeff)*IMU_pitch_p+IMU_filter_coeff*IMU_pitch;
knee=Knee_motor_position;
ankle=Ankle_motor_position;
dt=Iteration_time;
t=Iteration_time+t;


dknee=(1-filter_coeff)*dknee_p+filter_coeff*(knee-knee_pos_p)/dt;
dankle=(1-filter_coeff)*dankle_p+filter_coeff*(ankle-ankle_pos_p)/dt;


%%


%%



s_d=clamp((q_h_0-hip)/(q_h_0-q_h_min)*c,0,1);
s_a=clamp(1+(1-s_m)/(q_h_0-q_h_m)*(hip-q_h_0),0,1);

[knee_d_d,dknee_d_d,ankle_d_d,dankle_d_d]=winterFit(s_d*1000,a_knee,a_ankle,knee_ind,ankle_ind);
[knee_d_a,dknee_d_a,ankle_d_a,dankle_d_a]=winterFit(s_a*1000,a_knee,a_ankle,knee_ind,ankle_ind);

if state==1
    knee_d=knee_d_d;
    ankle_d=ankle_d_d;
    if s_d>=s_po
        state=2;
    end
elseif state==2
    knee_d=knee_d_d;
    ankle_d=ankle_d_d;
    if IMU_pitch_p<=hip
        state=3;
        q_h_m=hip;
        s_m=s_d;
    end
elseif state==3
    if IMU_pitch_p>hip
        knee_d=knee_p;
        ankle_d=ankle_p;
    else
        knee_d=knee_d_a;
        ankle_d=ankle_d_a;
    end
    if FC==0
        state=4;
    end
else
    knee_d=knee_d_a;
    ankle_d=ankle_d_a;
    if FC==1
        state=1;
    end
end




%%
if knee_ind_of_hip==1
    knee_des=knee_des_cte;
else
    knee_des=knee_d;
end
    

knee_des=clamp(knee_des,knee_pos_min_lim,knee_pos_max_lim);
knee_des=clamp(knee_des,dt*knee_vel_min_lim+knee_des_p,dt*knee_vel_max_lim+knee_des_p);
dknee_des=0;


if ankle_ind_of_hip==1
    ankle_des=ankle_des_cte;
else
    ankle_des=ankle_d;    
end

ankle_des=clamp(ankle_des,ankle_pos_min_lim,ankle_pos_max_lim);
ankle_des=clamp(ankle_des,dt*ankle_vel_min_lim+ankle_des_p,dt*ankle_vel_max_lim+ankle_des_p);
dankle_des=0;

Knee_torque_command=kp_knee*(knee_des-knee)+kd_knee*(dknee_des-dknee);
Ankle_torque_command=kp_ankle*(ankle_des-ankle)+kd_ankle*(dankle_des-dankle);
            





%%

%Storing persistent variables for next iteration
knee_pos_p=knee;
ankle_pos_p=ankle;
knee_des_p=knee_des;
ankle_des_p=ankle_des;
knee_p=knee_d;
ankle_p=ankle_d;
t_out=t;
time_p=time;
dknee_p=dknee;
dankle_p=dankle;
IMU_pitch_p=hip;
end %function

%%

function y=clamp(x,x1,x2)
y=min(x,max(x1,x2));
y=max(y,min(x1,x2));
end



function [knee_r,dknee_r,ankle_r,dankle_r]=winterFit(n,a_knee,a_ankle,knee_ind,ankle_ind)
m1=2;
m2=6;

c_n_knee=n-knee_ind;
j_knee=length(c_n_knee(c_n_knee>=0));
x=(n-knee_ind(j_knee))/(knee_ind(j_knee+1)-knee_ind(j_knee));
knee_r=a_knee(j_knee,1)*(x-1).^m1+a_knee(j_knee,2)*(x-1).^m2+a_knee(j_knee,4)*(x-1)+a_knee(j_knee,3);
dknee_r=(m1*a_knee(j_knee,1)*(x-1).^(m1-1)+m2*a_knee(j_knee,2)*(x-1).^(m2-1)+a_knee(j_knee,4))/(knee_ind(j_knee+1)-knee_ind(j_knee));


c_n_ankle=n-ankle_ind;
j_ankle=length(c_n_ankle(c_n_ankle>=0));
x=(n-ankle_ind(j_ankle))/(ankle_ind(j_ankle+1)-ankle_ind(j_ankle));
ankle_r=a_ankle(j_ankle,1)*(x-1).^m1+a_ankle(j_ankle,2)*(x-1).^m2+a_ankle(j_ankle,4)*(x-1)+a_ankle(j_ankle,3);
dankle_r=(m1*a_ankle(j_ankle,1)*(x-1).^(m1-1)+m2*a_ankle(j_ankle,2)*(x-1).^(m2-1)+a_ankle(j_ankle,4))/(ankle_ind(j_ankle+1)-ankle_ind(j_ankle));

end