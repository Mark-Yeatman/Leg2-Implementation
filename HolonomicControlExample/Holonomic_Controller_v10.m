function [Knee_torque_command, Ankle_torque_command, knee_des, ankle_des, ...
    t_out, hip, ankle_d_a, ankle_d_d, dt, state_out, s_a, s_d]...
    = Holonomic_Controller_v9(IMU_pitch, Knee_motor_position,...
    Knee_joint_position, Ankle_motor_position, Ankle_joint_position, Iteration,...
    Iteration_time, Knee_torque_sensor, Ankle_torque_sensor, Load_cell_x_force,...
    Load_cell_y_force, Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,...
    Load_cell_z_moment,kp_knee,kd_knee,kp_ankle,kd_ankle,ankle_ind_of_hip,...
    knee_ind_of_hip,ankle_des_cte,knee_des_cte,time,filter_coeff,IMU_filter_coeff,...
   q_h_0,q_h_min,c,FC,t_po,pushoff_init,s_a_sw,deltaT_init,s_td_fw,s_td_bw,kp_knee_td,...
   kd_knee_td,kp_ankle_td,kd_ankle_td,knee0,s_st_bw,n_td_bw,n_td_fw,t_td_bw,t_td_fw,...
   kp_knee_dual,kd_knee_dual,kp_ankle_dual,kd_ankle_dual,...
   ankle0_dual,ankle_m_dual,knee0_dual,ankle_m,n_r_2,n_sw)
%
% Inputs in addition to sensors:

% kp_ankle,kd_ankle,kp_knee,kd_knee are joint values

% ankle_ind_of_hip and knee_ind_of_hip are boolean values (true or false),
% --can be checkboxes in LabView with default TRUE.

% ankle_des_cte and knee_des_cte are desired angles for ankle and knee,
% controlled from LabView.


%% Variables Defined
%Persistent variables used to store data between iterations
persistent knee_pos_p ankle_pos_p ankle_des_p knee_des_p t time_p dknee_p dankle_p IMU_pitch_p q_h_m s_m state t0 ankle0_s knee0_s dknee0_s dankle0_s kp_knee0 kd_knee0 kp_ankle0 kd_ankle0

if isempty(time_p)
    time_p=time;
end
if isempty(knee_pos_p) || ((time-time_p)>1)
    knee_pos_p=0;
    ankle_pos_p=0;
    ankle_des_p=0;
    knee_des_p=10;
    knee_p=knee0;
    ankle_p=0;
    t=0;
    dknee_p=0;
    dankle_p=0;
    IMU_pitch_p=5;
    q_h_m=-10.99;
    s_m=c;
    state=1;
    t0=0;
    ankle0_s=0;
    knee0_s=knee_des_cte;
    dknee0_s=0;
    dankle0_s=0;
    kp_knee0=kp_knee;
    kd_knee0=kd_knee;
    kp_ankle0=kp_ankle;
    kd_ankle0=kd_ankle;
end


%% Winters Data Arrays
a_knee=[-21.9472500000000,4.24525000000000,21.6720000000000,0;20.9295000000000,-6.97650000000000,7.71900000000000,0;15.9362500000000,-0.998250000000000,18.6640000000000,25.8830000000000;-59.8412500000000,13.6422500000000,64.8630000000000,0;-19.4405000000000,-1.40250000000000,38.4100000000000,-47.2960000000000;45.1070000000000,-7.15300000000000,0.456000000000000,0;1.42525000000000,0.0657500000000000,2.21000000000000,3.24500000000000];
a_ankle=[5.50050000000000,-0.880500000000000,-4.60000000000000,0;-8.98500000000000,3.77020000000000,5.26600000000000,4.65120000000000;0.812500000000000,-0.475000000000000,7.74100000000000,2.81250000000000;-2.19975000000000,0.320750000000000,9.62000000000000,0;-14.7575000000000,0.878500000000001,-0.745000000000000,-24.2440000000000;23.0932000000000,-3.91420000000000,-19.9240000000000,0;-29.8020000000000,9.93400000000000,-0.0560000000000000,0;0.703500000000000,-0.234500000000000,-0.525000000000000,0;-2.60850000000000,0.869500000000000,1.21400000000000,0;-0.729000000000000,0.0190000000000000,0.580000000000000,-1.34400000000000];
knee_ind=[1 139 400 530 720 848 976 1001];
ankle_ind=[1 61 205 330 440 550 653 846 898 977 1001];

% n_sw=653;
% q_h_po=-8;
% n_r_2=433;


    
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
knee=Knee_joint_position;
ankle=Ankle_joint_position;
dt=Iteration_time;
t=Iteration_time+t;


dknee=(1-filter_coeff)*dknee_p+filter_coeff*(knee-knee_pos_p)/dt;
dankle=(1-filter_coeff)*dankle_p+filter_coeff*(ankle-ankle_pos_p)/dt;


%%


%%

kp_knee_f=kp_knee;
kd_knee_f=kd_knee;
kp_ankle_f=kp_ankle;
kd_ankle_f=kd_ankle;

s_d=clamp((q_h_0-hip)/(q_h_0-q_h_min)*c,0,1);
s_a=clamp(1+(1-s_m)/(q_h_0-q_h_m)*(hip-q_h_0),0,1);

[knee_d_d,dknee_d_d,ankle_d_d,dankle_d_d]=winterFit(s_d*1000,a_knee,a_ankle,knee_ind,ankle_ind);
[knee_d_a,dknee_d_a,ankle_d_a,dankle_d_a]=winterFit(s_a*1000,a_knee,a_ankle,knee_ind,ankle_ind);

if state==1 %stance/touchdown


%         knee_des=knee0;
%         dknee_des=0;
%         kp_knee_f=kp_knee_td;
%         kd_knee_f=kd_knee_td;
%         
%         ankle_des=0;
%         dankle_des=0;
%         kp_ankle_f=kp_ankle_td;
%         kd_ankle_f=kd_ankle_td;
        
        
        
        
    knee_d=knee0;
    ankle_d=0;

    [knee_des,temp]=spline3((t-t0)/deltaT_init,knee0_s,dknee0_s*deltaT_init,knee_d,0,0);
    dknee_d=temp/deltaT_init;
    dknee_des=0;  %no spline, since dknee_des is already zero in states 5 or 6 - even from state 4 (swing), zero velocity is best
    [ankle_des,temp]=spline3((t-t0)/deltaT_init,ankle0_s,dankle0_s*deltaT_init,ankle_d,0,0);
    dankle_d=temp/deltaT_init;
    dankle_des=0;
 
    
    kp_knee_f=spline3((t-t0)/deltaT_init,kp_knee0,0,kp_knee_td,0,0);
    kd_knee_f=spline3((t-t0)/deltaT_init,kd_knee0,0,kd_knee_td,0,0);
    kp_ankle_f=spline3((t-t0)/deltaT_init,kp_ankle0,0,kp_ankle_td,0,0);
    kd_ankle_f=spline3((t-t0)/deltaT_init,kd_ankle0,0,kd_ankle_td,0,0);
    
    if ankle>=ankle_m_dual
            state=2;
%             ankle0_2=(kp_ankle_dual-kp_ankle_td)/kp_ankle_dual*ankle;
    end
 
elseif state==2 %dual
    
    s_dual=(ankle-ankle_m_dual)/2;
        
    knee_des=spline3(s_dual,knee0,0,knee0_dual,0,0);
    dknee_des=0;
    
    kp_knee_f=spline3(s_dual,kp_knee_td,0,kp_knee_dual,0,0);
    kd_knee_f=spline3(s_dual,kd_knee_td,0,kd_knee_dual,0,0);
    
    kp_ankle_f=spline3(s_dual,kp_ankle_td,0,kp_ankle_dual,0,0);
    kd_ankle_f=spline3(s_dual,kd_ankle_td,0,kd_ankle_dual,0,0);
    
    ankle_des=spline3(s_dual,0,0,ankle0_dual,0,0);
    dankle_des=0;
    
     
    if s_d<0.2 || FC==0
        state=1;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=0;
        dankle0_s=0;
        kp_knee0=kp_knee_f;
        kd_knee0=kd_knee_f;
        kp_ankle0=kp_ankle_f;
        kd_ankle0=kd_ankle_f;
    end
    
    if ankle>=ankle_m
        state=3;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=0;
        dankle0_s=0;
        kp_knee0=kp_knee_f;
        kd_knee0=kd_knee_f;
        kp_ankle0=kp_ankle_f;
        kd_ankle0=kd_ankle_f;
    end

        
%     if FC==0 && s_a>0.99
%         state=5;
%         t0=t;
%         knee0_s=knee_des;
%         ankle0_s=ankle_des;
%         dknee0_s=0;
%         dankle0_s=0;
%     end
%     if FC==0 && s_d<0.01
%         state=6;
%         t0=t;
%         knee0_s=knee_des;
%         ankle0_s=ankle_des;
%         dknee0_s=0;
%         dankle0_s=0;
%     end
elseif state==3 %pushoff
    po_init=min(t_po,pushoff_init);
    n1=clamp(n_r_2+(n_sw-n_r_2)/t_po*(t-t0),n_r_2,n_sw);
    [knee_d,dknee_d,ankle_d,dankle_d]=winterFit(n1,a_knee,a_ankle,knee_ind,ankle_ind);
    dknee_d=dknee_d*(n_sw-n_r_2)/t_po;
    dankle_d=dankle_d*(n_sw-n_r_2)/t_po;
    
    [knee_des,temp]=spline3((t-t0)/po_init,knee0_s,dknee0_s*po_init,knee_d,dknee_d*po_init,1);
    dknee_des=temp/po_init;
    [ankle_des,temp]=spline3((t-t0)/po_init,ankle0_s,dankle0_s*po_init,ankle_d,dankle_d*po_init,1);
    dankle_des=temp/po_init;
    
    kp_knee_f=spline3((t-t0)/deltaT_init,kp_knee0,0,kp_knee,0,0);
    kd_knee_f=spline3((t-t0)/deltaT_init,kd_knee0,0,kd_knee,0,0);
    kp_ankle_f=spline3((t-t0)/deltaT_init,kp_ankle0,0,kp_ankle,0,0);
    kd_ankle_f=spline3((t-t0)/deltaT_init,kd_ankle0,0,kd_ankle,0,0);
 
    if s_a>s_a_sw && t>=(t0+t_po)
        state=4;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        t0=t;
        dknee0_s=dknee_des;
        dankle0_s=dankle_des;
        kp_knee0=kp_knee;
        kd_knee0=kd_knee;
        kp_ankle0=kp_ankle;
        kd_ankle0=kd_ankle;
    end
elseif state==4 %swing
 
    knee_d=knee_d_a;
    ankle_d=ankle_d_a;

    [knee_des,dknee_d]=spline3((t-t0)/deltaT_init,knee0_s,dknee0_s*deltaT_init,knee_d,dknee_d_a*deltaT_init,1);
    [dknee_des,temp]=spline3((t-t0)/deltaT_init,dknee_d/deltaT_init,0,0,0,0);
    [ankle_des,dankle_d]=spline3((t-t0)/deltaT_init,ankle0_s,dankle0_s*deltaT_init,ankle_d,dankle_d_a*deltaT_init,1);
    [dankle_des,temp]=spline3((t-t0)/deltaT_init,dankle_d/deltaT_init,0,0,0,0);
    
    kp_knee_f=spline3((t-t0)/deltaT_init,kp_knee0,0,kp_knee,0,0);
    kd_knee_f=spline3((t-t0)/deltaT_init,kd_knee0,0,kd_knee,0,0);
    kp_ankle_f=spline3((t-t0)/deltaT_init,kp_ankle0,0,kp_ankle,0,0);
    kd_ankle_f=spline3((t-t0)/deltaT_init,kd_ankle0,0,kd_ankle,0,0);
    
    if FC==0 && s_a>s_td_fw
        state=5;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_d/deltaT_init;
        dankle0_s=dankle_d/deltaT_init;
    end
   
    if FC==0 && s_d>s_td_bw
        state=6;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_d/deltaT_init;
        dankle0_s=dankle_d/deltaT_init;
    end
    
    if FC==1 && knee<30 % && s_a>s_a_sw
        state=1;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_d/deltaT_init;
        dankle0_s=dankle_d/deltaT_init;
        kp_knee0=kp_knee;
        kd_knee0=kd_knee;
        kp_ankle0=kp_ankle;
        kd_ankle0=kd_ankle;
    end
%     
elseif state==5 %touchdown-forward
%     n_td_fw=870;
%     t_td_fw=0.07;
    deltaT_td=t_td_fw/2;
    
    n1=clamp(n_td_fw+(1000-n_td_fw)/t_td_fw*(t-t0),0,1000);
    [knee_d,dknee_d,ankle_d,dankle_d]=winterFit(n1,a_knee,a_ankle,knee_ind,ankle_ind);
    dknee_d=dknee_d*(1000-n_td_fw)/t_td_fw;
    dankle_d=dankle_d*(1000-n_td_fw)/t_td_fw;
    
    kp_knee_f=spline3((t-t0)/deltaT_td,kp_knee,0,kp_knee_td,0,0);
    kd_knee_f=spline3((t-t0)/deltaT_td,kd_knee,0,kd_knee_td,0,0);
    kp_ankle_f=spline3((t-t0)/deltaT_td,kp_ankle,0,kp_ankle_td,0,0);
    kd_ankle_f=spline3((t-t0)/deltaT_td,kd_ankle,0,kd_ankle_td,0,0);
    
    [knee_des,temp]=spline3((t-t0)/deltaT_td,knee0_s,dknee0_s*deltaT_td,knee_d,dknee_d*deltaT_td,1);
    dknee_des=0;
    [ankle_des,temp]=spline3((t-t0)/deltaT_td,ankle0_s,dankle0_s*deltaT_td,ankle_d,dankle_d*deltaT_td,1);
    dankle_des=0;
    
    knee_des=max(knee_des,knee0);

            

    if FC==1 && s_d>0.1 && t>=(t0+t_td_fw)
        state=1;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_des;
        dankle0_s=dankle_des;
        kp_knee0=kp_knee_f;
        kd_knee0=kd_knee_f;
        kp_ankle0=kp_ankle_f;
        kd_ankle0=kd_ankle_f;
    end
    if FC==0 && s_d>0.1
        state=4;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_des;
        dankle0_s=dankle_des;
        kp_knee0=kp_knee_f;
        kd_knee0=kd_knee_f;
        kp_ankle0=kp_ankle_f;
        kd_ankle0=kd_ankle_f;
    end
    
else %touchdown-backward
%     n_td_bw=870;
%     t_td_bw=0.15;
    deltaT_td=t_td_bw/2;
    
    n1=clamp(n_td_bw+(1000-n_td_bw)/t_td_bw*(t-t0),0,1000);
    [knee_d,dknee_d,ankle_d,dankle_d]=winterFit(n1,a_knee,a_ankle,knee_ind,ankle_ind);
    dknee_d=dknee_d*(1000-n_td_bw)/t_td_bw;
    dankle_d=dankle_d*(1000-n_td_bw)/t_td_bw;
    
    kp_knee_f=spline3((t-t0)/deltaT_td,kp_knee,0,kp_knee_td,0,0);
    kd_knee_f=spline3((t-t0)/deltaT_td,kd_knee,0,kd_knee_td,0,0);
    kp_ankle_f=spline3((t-t0)/deltaT_td,kp_ankle,0,kp_ankle_td,0,0);
    kd_ankle_f=spline3((t-t0)/deltaT_td,kd_ankle,0,kd_ankle_td,0,0);
    
    [knee_des,temp]=spline3((t-t0)/deltaT_td,knee0_s,dknee0_s*deltaT_td,knee_d,dknee_d*deltaT_td,1);
    dknee_des=0;
    [ankle_des,temp]=spline3((t-t0)/deltaT_td,ankle0_s,dankle0_s*deltaT_td,ankle_d,dankle_d*deltaT_td,1);
    dankle_des=0;
    
    knee_des=max(knee_des,knee0);
           

    if FC==1 && s_d<s_st_bw && t>=(t0+t_td_bw)
        state=1;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_des;
        dankle0_s=dankle_des;
        kp_knee0=kp_knee_f;
        kd_knee0=kd_knee_f;
        kp_ankle0=kp_ankle_f;
        kd_ankle0=kd_ankle_f;
    end
    if FC==0 && s_d<s_st_bw
        state=4;
        t0=t;
        knee0_s=knee_des;
        ankle0_s=ankle_des;
        dknee0_s=dknee_des;
        dankle0_s=dankle_des;
        kp_knee0=kp_knee_f;
        kd_knee0=kd_knee_f;
        kp_ankle0=kp_ankle_f;
        kd_ankle0=kd_ankle_f;
    end
end




%%
if knee_ind_of_hip==1
    knee_des=knee_des_cte;
    dknee_des=0;
end
    

knee_des=clamp(knee_des,knee_pos_min_lim,knee_pos_max_lim);
knee_des=clamp(knee_des,dt*knee_vel_min_lim+knee_des_p,dt*knee_vel_max_lim+knee_des_p);
dknee_des=clamp(dknee_des,knee_vel_min_lim,knee_vel_max_lim);

if ankle_ind_of_hip==1
    ankle_des=ankle_des_cte;
    dankle_des=0;
end

ankle_des=clamp(ankle_des,ankle_pos_min_lim,ankle_pos_max_lim);
ankle_des=clamp(ankle_des,dt*ankle_vel_min_lim+ankle_des_p,dt*ankle_vel_max_lim+ankle_des_p);
dankle_des=clamp(dankle_des,ankle_vel_min_lim,ankle_vel_max_lim);


Knee_torque_command=kp_knee_f*(knee_des-knee)+kd_knee_f*(dknee_des-dknee);
Ankle_torque_command=kp_ankle_f*(ankle_des-ankle)+kd_ankle_f*(dankle_des-dankle);
            





%%

%Storing persistent variables for next iteration
knee_pos_p=knee;
ankle_pos_p=ankle;
knee_des_p=knee_des;
ankle_des_p=ankle_des;
% knee_p=knee_d;
% ankle_p=ankle_d;

time_p=time;
dknee_p=dknee;
dankle_p=dankle;
IMU_pitch_p=hip;

state_out=state;
t_out=t;
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


function [y,dy]=spline3(x,y0,dy0,y1,dy1,sat_spd)
%x between 0 and 1

a0=y0;
a1=dy0;
a2=3*y1-3*y0-2*dy0-dy1;
a3=-2*y1+2*y0+dy0+dy1;

if x<0
    y=y0;
    if sat_spd
        dy=dy0;
    else
        dy=0;
    end
elseif x>1
    y=y1;
    if sat_spd
        dy=dy1;
    else
        dy=0;
    end
else
    y=a0+a1*x+a2*x.^2+a3*x.^3;
    dy=a1+2*a2*x+3*a3*x.^2;
end

end
