function [Knee_torque_command, Ankle_torque_command, deltaL, hip_pos, Esys, Esys_integrate_out,...
        U_S_KNEE, U_S_ANKLE, COPFX, U_LIN_DAMP_A, U_STOP_K, U_STOP_A, U_PBC_K,...
        U_PBC_A, knee_des_out,ankle_des_out, foot_contact, stance, swing,phase_var_out,...
        IMU_LIVE_OUT,StanceGain,SwingGain,knee_joint_vel,ankle_joint_vel,hip_vel,PushOffOut]...
= SLIP_KPBC(IMU_pitch, Knee_joint_position, Ankle_joint_position, ...
            time_in, dt, ...
            Load_cell_x_force, Load_cell_y_force, Load_cell_z_force,...
            Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment,...
            kp_knee, kd_knee, kp_ankle, kd_ankle, ankle_des_in, knee_des_in,...
            vel_filter_coeff, KPBC_filter_coeff,...
            SLIP_ON, lt, k, d, L0, ...
            KPBC_ON, KPBC_max_torque_rate , pbc_gain_knee, md, Eref,...
            knee_stop_low, knee_stop_high, ankle_stop_low, ankle_stop_high,...
            max_torque,...
            v0, Fric_Comp, F_thresh,...
            Command_State,...
            KPBC_max_torque, Joint_Bio_Sat)

%Inputs in addition to sensors:

% Knee_joint_position, Ankle_joint_position are in degrees, in a
% "biomechanics frame"

% t, dt, t_out are in seconds

    %% Persistent Variable Definitions
%Persistent variables used to store data between iterations
persistent hip_pos_prev...
           knee_pos_prev...
           ankle_pos_prev...
           knee_vel_prev...
           ankle_vel_prev...
           hip_vel_prev...
           Esys_integrate...
           t_switched_phase...
           Stance...
           Swing...
           ForceCount...
           Hip_Pos_Max...
           Hip_Pos_Min...
           IMU_LIVE...
           PushOff...
           u_pbc_knee_prev...
           u_pbc_ankle_prev...
           %t_switched_swing... %implement swing setpoint switching hold?
       
    if isempty(knee_pos_prev)
        hip_pos_prev = IMU_pitch;
        knee_pos_prev = 0;
        ankle_pos_prev = 0;
        knee_vel_prev = 0;
        ankle_vel_prev = 0;
        hip_vel_prev = 0;
        IMU_pitch_prev = 0;
        Esys_integrate  = 0;
        u_pbc_knee_prev = 0;
        u_pbc_ankle_prev = 0;
        t_switched_phase = time_in;
        %t_switched_swing = time_in;
        Stance = true;
        Swing = false;
        ForceCount = 0;
        Edis = 0;
        Hip_Pos_Max  = 23;
        Hip_Pos_Min = -20;
        IMU_LIVE = true;
        PushOff = false;
    end

    %% Initialization and Hard Coded Values
    Ms = single(6.991429/2.205); % converted to kg
    Mf = single(0.55 + 1.643130/2.205); %carbon fiber foot + mechanism, kg
    %lt = single(0.3733); %meters
    ls = single(0.3733); %meters
    la = single(0.0628); %meters
    lfx = single(-0.0071882); %meters
    lfy = single(0.00338582); %meters
    lc = single(2.00914e-5); %meters
    g = 9.81; %meters/s^2;
    

    
    COPFX = single(0);
    U_LIN_DAMP_A = single(0);
    U_S_KNEE = single(0);
    U_S_ANKLE = single(0);
    U_PBC_K = single(0);
    U_PBC_A = single(0);
    
    % Position/velocity hard limits
    ANKLE_POS_MIN_LIM  = -35;
    ANKLE_POS_MAX_LIM = 35;
    ANKLE_VEL_MIN_LIM = -200; %deg/s
    ANKLE_VEL_MAX_LIM = 200; %deg/s
    KNEE_POS_MIN_LIM = 2;
    KNEE_POS_MAX_LIM = 105;
    KNEE_VEL_MIN_LIM = -400; %deg/s
    KNEE_VEL_MAX_LIM = 400; %deg/s
    
    %Winters Data Arrays
    a_knee=[-21.9472500000000,4.24525000000000,21.6720000000000,0;20.9295000000000,-6.97650000000000,7.71900000000000,0;15.9362500000000,-0.998250000000000,18.6640000000000,25.8830000000000;-59.8412500000000,13.6422500000000,64.8630000000000,0;-19.4405000000000,-1.40250000000000,38.4100000000000,-47.2960000000000;45.1070000000000,-7.15300000000000,0.456000000000000,0;1.42525000000000,0.0657500000000000,2.21000000000000,3.24500000000000];
    a_ankle=[5.50050000000000,-0.880500000000000,-4.60000000000000,0;-8.98500000000000,3.77020000000000,5.26600000000000,4.65120000000000;0.812500000000000,-0.475000000000000,7.74100000000000,2.81250000000000;-2.19975000000000,0.320750000000000,9.62000000000000,0;-14.7575000000000,0.878500000000001,-0.745000000000000,-24.2440000000000;23.0932000000000,-3.91420000000000,-19.9240000000000,0;-29.8020000000000,9.93400000000000,-0.0560000000000000,0;0.703500000000000,-0.234500000000000,-0.525000000000000,0;-2.60850000000000,0.869500000000000,1.21400000000000,0;-0.729000000000000,0.0190000000000000,0.580000000000000,-1.34400000000000];
    knee_ind=[1 139 400 530 720 848 976 1001];
    ankle_ind=[1 61 205 330 440 550 653 846 898 977 1001];
    
    %% State Calculations
    %Assign joint positions
    hip_pos = IMU_pitch;
    knee_pos =  Knee_joint_position; 
    ankle_pos = Ankle_joint_position;
    
    %Calculate joint velocities usindg exponential smoothing filter
    %https://en.wikipedia.org/wiki/Exponential_smoothing
    hip_vel = (1-vel_filter_coeff)*hip_vel_prev + vel_filter_coeff*(hip_pos-hip_pos_prev)/dt;
    knee_vel = (1-vel_filter_coeff)*knee_vel_prev + vel_filter_coeff*(knee_pos-knee_pos_prev)/dt;
    ankle_vel = (1-vel_filter_coeff)*ankle_vel_prev + vel_filter_coeff*(ankle_pos-ankle_pos_prev)/dt;
    
    %% Determine foot contact
    if norm([Load_cell_x_force, Load_cell_y_force, Load_cell_z_force],2) > F_thresh
        ForceCount = ForceCount+1;
        if ForceCount > 2
            FootContact = true;
        else
            FootContact = false;
        end
    else
        FootContact = false;
        ForceCount = 0;   
    end
             
    StanceGain = 0;
    SwingGain = 0;
    t_hold = 0.3;
    
    %Human Leg as a Robot states
    x = zeros(4,1); 
    x(1) =  knee_pos; %MAKE SURE TO MANAGE BIOMECHANICS VS RIGHT-HAND-RULE AXIS ORIENTATION
    x(2) =  ankle_pos;
    x(3) =  knee_vel; 
    x(4) =  ankle_vel;
       
    %Notes:
    % Knee axis treated as origin.
    % Orientation of knee and ankle rotation is biomechanical

    %variables:
    % x - 4x1 array [knee pos, ankle pos, knee vel, ankle vel]
    % lt -thigh length 
    % ls -shank length
    % lfx, lfy - foot CoM x,y
    % lc - load cell y dist from ankle axis
    % la - ankle axis y dist to bottom of foot
    % COPfx - COP distance along foot from ankle axis projection, to be
    %   computed from load cell
    % Ms - shank mass
    % Mf - foot mass

    % COPfx calc
    n = [0,1,0]';
    MA = [Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment]';
    F = [Load_cell_x_force, Load_cell_y_force, Load_cell_z_force]';
    WC = [F;MA];   
    OA = la*n;
    M0 = (cross(-OA,F)+MA);
    OC = cross(n,M0)/dot(F,n);
    COPfx = OC(1);        
    if ~FootContact
        %(-2,8) inches
        COPfx = (0.2032 - 0.0508)/2; %meters
    end
    COPfx = Saturate(COPfx,-0.0508,0.2032); %meters
    COPFX = COPfx;
    
    params = [Ms,Mf,Isz,Ifz,lt,ls,lc,la,rs,rfy,rfx,COPfx,g];
   
    M = M_func(x,params);
    C = C_func(x,params);
    G = G_func(x,params);    
    L = L_func(x,params);
    J_L = J_L_func(x,params);
    J_dot_L = J_dot_L_func(x,params);    
    Theta = Theta_func(x,params);
    J_Theta = J_Theta_func(x,params);
    J_dot_Theta = J_dot_Theta_func(x,params);   
    J_z = [J_L;J_Theta];
    J_z_dot = [J_dot_L;J_dot_Theta];      
    J_C = J_C_func(x,params);   
    
    qdot = x(3:4);
    deltaL = L - L0;
    z_dot = J_z*qdot;
    Ldot = z_dot(1);
    if md<10
       md = 10;
    elseif md>10000
        md = 10000;
    end
    M_z = diag([md,md]);
    
    %Calculate system energy
    Espring = 1/2*k*deltaL^2; %virtual spring potential energy
    KE = 1/2*md*(Ldot^2);     %virtual mass kinetic energy
    Esys = Espring + KE;
    
    %Phase Variable 
    %s_a=clamp(1+(1-s_m)/(q_h_0-q_h_m)*(hip-q_h_0),0,1);
    phase_var_out =  0.6 + 0.4*(hip_pos-Hip_Pos_Min)/(Hip_Pos_Max-Hip_Pos_Min);
    phase_var_out = clamp(phase_var_out,0.6,1);
    
    %IMU State check
    IMU_LIVE = abs(hip_pos_prev-hip_pos)<1;
    
    %% Phase State management, set autodetection or command stance/swing
    if int8(Command_State) == int8(0)
        %Determine stance vs swing, only allow phase state switching every t_hold seconds
        if Stance && ~FootContact %The state needs to be switched
            if t_hold < (time_in - t_switched_phase) %Check the last time we switched, if its too fast, wait
               %Switch to SWING
               t_switched_phase = time_in;
               Swing = true;
               Stance = false;
               %reset the energy intergration variables
               %Esys_integrate = 0;
               %Edis = 0; 
               u_pbc_knee_prev = 0;
               u_pbc_ankle_prev = 0;
            end
        elseif Swing && FootContact %ditto, but for the other case
            if t_hold < (time_in - t_switched_phase) 
                %Switch to STANCE
                t_switched_phase = time_in;
                Swing = false;
                Stance = true;
            end
        end
        %Manage torque smoothing when switching phases
        if Stance
            t_norm = (time_in - t_switched_phase)/(t_hold*0.8);
            StanceGain = min(t_norm,1);
            SwingGain = 1-StanceGain; 
        elseif Swing
            t_norm = (time_in - t_switched_phase)/(t_hold*0.8);
            SwingGain = min(t_norm,1);
            StanceGain = 1-SwingGain; 
        end
    elseif int8(Command_State) == int8(1) %Stance
        Stance = true;
        Swing = false;
        Hip_Pos_Max = 23;
        Hip_Pos_Min = - 20;
        StanceGain = 1;
        SwingGain = 0;
    elseif int8(Command_State) == int8(2) %Swing
        Stance = false;
        Swing = true;
        Hip_Pos_Max = 23;
        Hip_Pos_Min = - 20;
        StanceGain = 0;
        SwingGain = 1;
    end

    %% Torque computation
    if SLIP_ON       
        Knee_torque_command_stance = 0;
        Ankle_torque_command_stance = 0;
        %% Stance                
        %Virtual spring
        z_tilde = [Theta-Theta_ref;L-L0];
        K = diag([k,k]);
        D = diag([d,d]);
        M_J_L_inv = pinv((J_L/M));
        u_s = C*qdot + G + M_J_L_inv*(J_z_dot*qdot + M_z\(-K*z_tilde - D*z_dot));

        %MAKE SURE TO MANAGE BIOMECHANICS VS RIGHT-HAND-RULE AXIS ORIENTATION
        U_S_KNEE = u_s(1);
        U_S_ANKLE = u_s(2);

        Knee_torque_command_stance  = Knee_torque_command_stance + u_s(1); 
        Ankle_torque_command_stance = Ankle_torque_command_stance + u_s(2);

        if KPBC_ON %Also use energy tracking controller   
            Omega = diag([1,0]);
            kappa = pbc_gain_knee;
            u_r  = -kappa*(Esys - Eref)*M_J_L_inv*Omega*Ldot; 

            U_PBC_K = u_r(1); 
            U_PBC_A = u_r(2);
                         
            %Moving average filter and saturation laws on value and rate of torque       
            %Absolute saturation
            if KPBC_max_torque>0
                U_PBC_K = Saturate(U_PBC_K,-KPBC_max_torque,KPBC_max_torque);
                U_PBC_A = Saturate(U_PBC_A,-KPBC_max_torque,KPBC_max_torque);
            end           
            %Rate saturation
            if KPBC_max_torque_rate>0
                dk = (U_PBC_K-u_pbc_knee_prev)/dt;
                da = (U_PBC_A-u_pbc_ankle_prev)/dt;
                dk = Saturate(dk, -KPBC_max_torque_rate, KPBC_max_torque_rate);
                da = Saturate(da, -KPBC_max_torque_rate, KPBC_max_torque_rate);
                U_PBC_K = dk*dt+u_pbc_knee_prev;
                U_PBC_A = da*dt+u_pbc_ankle_prev;
            end           
            %Exponential smoothing of torque
            U_PBC_K = (1-KPBC_filter_coeff)*u_pbc_knee_prev + KPBC_filter_coeff*U_PBC_K;
            U_PBC_A = (1-KPBC_filter_coeff)*u_pbc_ankle_prev + KPBC_filter_coeff*U_PBC_A;
            
            %Biomemetic power saturation
            if Joint_Bio_Sat               
                Pow_Knee = knee_vel*U_PBC_K ;
                Pow_Ankle = ankle_vel*U_PBC_A;
                if sign(Pow_Ankle) == -1
                    U_PBC_A = 0;
                end
                if sign(Pow_Knee) == 1
                    U_PBC_K = 0;
                end
            end
                      
        end
        %% Swing
        if IMU_LIVE
            %Use hip pos as phase variable to index winters data
            [knee_des_out ,~,ankle_des_out,~]=winterFit(phase_var_out*1000,a_knee,a_ankle,knee_ind,ankle_ind);

            %follow the trjactory via PD controller
            [Knee_torque_command_swing,Ankle_torque_command_swing] = PDControl(...
            kp_knee, kd_knee, knee_des_out, knee_pos, knee_vel, ...
            kp_ankle, kd_ankle, ankle_des_out, ankle_pos, ankle_vel);
       
        else
            knee_des_out = knee_pos;
            ankle_des_out = ankle_pos;
            [Knee_torque_command_swing,Ankle_torque_command_swing] = PDControl(...
            kp_knee, kd_knee, knee_des_out, knee_pos, knee_vel, ...
            kp_ankle, kd_ankle, ankle_des_out, ankle_pos, ankle_vel);
        end
        
        %Smooth between stance and swing torques
        Knee_torque_command = StanceGain*Knee_torque_command_stance + SwingGain*Knee_torque_command_swing;
        Ankle_torque_command = StanceGain*Ankle_torque_command_stance + SwingGain*Ankle_torque_command_swing;
        
    else
        IMU_LIVE = true;
        %Use Setpoint PD control
        knee_des_out = knee_des_in;
        ankle_des_out = ankle_des_in;
        [Knee_torque_command,Ankle_torque_command] = PDControl(...
            kp_knee, kd_knee, knee_des_out, knee_pos, knee_vel, ...
            kp_ankle, kd_ankle, ankle_des_out, ankle_pos, ankle_vel);        
    end
    
    %% Friction Compensation
    %v0 = 0.01;
    tau0 = 5;
    if Fric_Comp
        Knee_torque_command = friction_compensation(knee_vel,Knee_torque_command ,v0,tau0);
        Ankle_torque_command = friction_compensation(ankle_vel,Ankle_torque_command ,v0,tau0);
    end
    
    %% Virtual hard stops
    u_stop = zeros(2,1);
    if knee_stop_low < knee_stop_high && ankle_stop_low < ankle_stop_high
        knee_limits = [knee_stop_low, knee_stop_high];  %deg, need to make sure its ordered
        ankle_limits = [ankle_stop_low, ankle_stop_high];
    else
        knee_limits = single([KNEE_POS_MIN_LIM, KNEE_POS_MAX_LIM]);  %deg, need to make sure its ordered
        ankle_limits = single([ANKLE_POS_MIN_LIM, ANKLE_POS_MAX_LIM]);
    end

    if knee_pos < knee_limits(1)
        u_stop(1) = -kp_knee*(knee_pos - knee_limits(1)) - kd_knee*(knee_vel);
    end
    if knee_limits(2) < knee_pos 
        u_stop(1) = -kp_knee*(knee_pos - knee_limits(2)) - kd_knee*(knee_vel);
    end
    if ankle_pos < ankle_limits(1)
        u_stop(2) = -kp_ankle*(ankle_pos - ankle_limits(1)) - kd_ankle*(ankle_vel);
    end
    if ankle_limits(2) < ankle_pos
        u_stop(2) = -kp_ankle*(ankle_pos - ankle_limits(2)) - kd_ankle*(ankle_vel);
    end
    U_STOP_K = u_stop(1);
    U_STOP_A = u_stop(2);
    Knee_torque_command = Knee_torque_command + u_stop(1);
    Ankle_torque_command = Ankle_torque_command + u_stop(2);
    
    %% Saturate torque output
    if max_torque>0
        Knee_torque_command = Saturate(Knee_torque_command,-max_torque,max_torque);
        Ankle_torque_command = Saturate(Ankle_torque_command,-max_torque,max_torque);
    end
        
    %% Persistent/output variable management
    Esys_integrate_out = Esys_integrate;
    knee_pos_prev = knee_pos;
    ankle_pos_prev = ankle_pos;
    knee_vel_prev=knee_vel;
    ankle_vel_prev=ankle_vel;
    IMU_pitch_prev=IMU_pitch;
    hip_vel_prev = hip_vel;
    hip_pos_prev = hip_pos;
    u_pbc_knee_prev = U_PBC_K;
    u_pbc_ankle_prev = U_PBC_A;
    knee_joint_vel = single(knee_vel);
    ankle_joint_vel = single(ankle_vel);
    foot_contact = single(FootContact);
    stance = single(Stance);
    swing = single(Swing);
    IMU_LIVE_OUT = single(IMU_LIVE);
    PushOffOut = single(PushOff);
end

%% Helper functions
function [Knee_torque_command,Ankle_torque_command] = PDControl(kp_k, kd_k, q_kstar, q_k, qdot_k, kp_a, kd_a, q_astar, q_a, qdot_a)
    %This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero.
    Knee_torque_command = kp_k*(q_kstar - q_k) + kd_k*(-qdot_k);
    Ankle_torque_command = kp_a*(q_astar - q_a) + kd_a*(-qdot_a);
end
function y = Saturate(x,L1,L2)
    %Function to prevent the desired joint angles from changing to fast. 
    %Works via saturation
    y = min(x,max(L1,L2));
    y = max(y,min(L1,L2));
end
function tau_com_w_fric=friction_compensation(v,tau_com,v0,tau0)

f_k_n=-1.7;
f_s_n=-2.8;

f_s_p=1.6;
f_k_p=0.8;

r_v_p=sat_one(v,0,v0);
r_v_n=sat_one(v,-v0,0);
r_t_p=sat_one(tau_com,0,tau0);
r_t_n=sat_one(tau_com,-tau0,0);

tau_fric=r_v_p*f_k_p+(1-r_v_n)*f_k_n+(1-r_v_p)*r_v_n*(r_t_p*f_s_p+(1-r_t_n)*f_s_n);

tau_com_w_fric=tau_com+tau_fric;
end
function y=sat_one(x,a,b)
if x<a
    y=0;
elseif x<b
    y=(x-a)/(b-a);
else
    y=1;
end
end
function [knee_r,dknee_r,ankle_r,dankle_r] = winterFit(n,a_knee,a_ankle,knee_ind,ankle_ind)
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
function y = clamp(x,x1,x2)
y=min(x,max(x1,x2));
y=max(y,min(x1,x2));
end