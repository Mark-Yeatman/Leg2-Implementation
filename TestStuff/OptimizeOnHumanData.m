set(0,'DefaultFigureWindowStyle','docked')

%Kyles Data is on my desktop at -> "C:\Users\mxy110230\Documents\KylesData\InclineExperiment.mat"

Subjects = fieldnames(Gaitcycle);

trial = 's1i0';

for j=1:length(Subjects)
    
    s = Subjects{j};
    
    lf = Limbs.(s).FromMarkers.BallOFoot/1000; %magic numbers from from page 83 in winter's biomechanics
    la = Gaitcycle.(s).subjectdetails{3,2}/1000*(0.039)*1.2; 
    ls = Limbs.(s).FromMarkers.Shank/1000;
    lt = Limbs.(s).FromMarkers.Thigh/1000; 
    
    m_0 = Gaitcycle.AB06.subjectdetails{4,2};
    k_0 = 15250-300; 
    %k = k_0 * Gaitcycle.(subject).subjectdetails{4,2} / m_0 ;
    d = 0;
    fc = lf;

    links = [lt,ls,la,lf];
    
    ankle_trajectory = Gaitcycle.(s).(trial).kinematics.jointangles.right.ankle.x_mean;
    knee_trajectory = Gaitcycle.(s).(trial).kinematics.jointangles.right.knee.x_mean;
    hip_trajectory = Gaitcycle.(s).(trial).kinematics.jointangles.right.hip.x_mean;

    x = [hip_trajectory, knee_trajectory, ankle_trajectory];
    
    xstart = x(1,:);
    xstart(2) = -xstart(2);
    L0 = (fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cosd(xstart(2))+la.*ls.*cosd(xstart(3))+la.*lt.*cosd(xstart(2)+xstart(3))+(-2).*fc.*ls.*sind(xstart(3))+(-2).*fc.*lt.*sind(xstart(2)+xstart(3))).^(1/2);
    
    ti = 1:length(ankle_trajectory);
    hipspline = spline(ti, hip_trajectory);
    p_der=fnder(hipspline,1);
    hip_vel = ppval(p_der,ti);
    
    kneespline = spline(ti, knee_trajectory);
    p_der=fnder(kneespline,1);
    knee_vel = ppval(p_der,ti);
    
    anklespline = spline(ti, ankle_trajectory);
    p_der=fnder(anklespline,1);
    ankle_vel = ppval(p_der,ti);
    
    cut_v = 0.6;
    cut_i = cut_v*length(ankle_trajectory);
    
    scale = Gaitcycle.(s).subjectdetails{4,2}/1000; %turns N*mm/kg to N*m
        
    steps = size(Gaitcycle.(s).(trial).kinematics.jointangles.right.ankle.x,2);
    
    tau_k_data = zeros(steps*cut_i,1);%Gaitcycle.(s).(trial).kinetics.jointmoment.right.knee.x_mean(1:cut_i)*scale;
    tau_a_data = zeros(steps*cut_i,1);%Gaitcycle.(s).(trial).kinetics.jointmoment.right.ankle.x_mean(1:cut_i)*scale;
    
    count = 1;
    for k = 1:steps
        for i = 1:cut_i
            tau_a_data(count) = Gaitcycle.(s).(trial).kinetics.jointmoment.right.ankle.x(i,k);
            tau_k_data(count) = Gaitcycle.(s).(trial).kinetics.jointmoment.right.knee.x(i,k);
    
            count = count+1;
        end
    end

    count = 1;
    clear kcvx
    cvx_begin 
        variables kcvx kL0cvx
        variables tau_k_slip(steps*cut_i) tau_a_slip(steps*cut_i)
        minimize 10*norm(tau_a_slip + tau_a_data) + 1*norm(tau_k_slip - tau_k_data)
        subject to         
            0 <= kcvx;
            for k = 1:steps
                for i = 1:cut_i
                    ankle = Gaitcycle.(s).(trial).kinematics.jointangles.right.ankle.x(i,k);
                    knee = Gaitcycle.(s).(trial).kinematics.jointangles.right.knee.x(i,k);
                    hip = Gaitcycle.(s).(trial).kinematics.jointangles.right.hip.x(i,k);
                    [tau_k_slip(count), tau_a_slip(count)] == SLIP_KPBC_forsim(...
                    hip, knee, ankle, 0, 0, 0, lf, la, ls, lt, kcvx, 0, L0, 0, 0, kL0cvx);
                                      
                    count = count+1;
                end
            end
    cvx_end
    
    K(j) = kcvx;  
    l0(j) = kL0cvx/kcvx;
    %Knee_shift(j) = tau_k_shift;
    %Ankle_shift(j) = tau_a_shift;
    
    figure('Name',s,'NumberTitle','off')
    subplot(2,2,1)
    ilist = linspace(0,cut_v,cut_i);
    plot(tau_k_data)
    hold on
    plot(tau_k_slip)
    title('Knee Torque Comparison')
    xlabel('Normalized Time')
    ylabel('Torque(N m)')
    legend('Human','SLIP')

    subplot(2,2,2)
    plot(tau_a_data)
    hold on 
    plot(-tau_a_slip)
    title('Ankle Torque Comparison')
    xlabel('Normalized Time')
    ylabel('Torque(N m)')
    legend('Human','SLIP')
    
end

figure('Name','K','NumberTitle','off')
ilist = linspace(0,cut_v,cut_i);
scatter(1:length(Subjects), K)
title('Knee Torque Comparison')
xlabel('Subject Number')
ylabel('K')


set(0,'DefaultFigureWindowStyle','normal')


%animate(@drawProsthesis, x,[1:length(ankle_trajectory)],links)
