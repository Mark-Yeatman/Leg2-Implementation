set(0,'DefaultFigureWindowStyle','docked')

%Kyles Data is on my desktop at -> "C:\Users\mxy110230\Documents\KylesData\InclineExperiment.mat"

Subjects = fieldnames(Gaitcycle); %{'AB06'};%
trial = 's1i0';

for j=1:length(Subjects)
    subject = Subjects{j};
    lf = 0.6*Gaitcycle.(subject).subjectdetails{3,2}/1000 *(0.152); %magic numbers from from page 83 in winter's biomechanics
    la = Gaitcycle.(subject).subjectdetails{3,2}/1000 *(0.039); 
    ls = Gaitcycle.(subject).subjectdetails{3,2}/1000 *(0.285)-la; 
    lt = Gaitcycle.(subject).subjectdetails{3,2}/1000 *(0.530)-ls-la; 
    
    m_0 = Gaitcycle.AB06.subjectdetails{4,2};
    k_0 = 15250; 
    k = k_0;% K(j);
    d = 0;%D(j);
    fc = lf;

    links = [lt,ls,la,lf];
    %Calls the control function. The control function houses all of the control
    %logic and mathematics.

    ankle_trajectory = Gaitcycle.(subject).(trial).kinematics.jointangles.right.ankle.x_mean;
    knee_trajectory = Gaitcycle.(subject).(trial).kinematics.jointangles.right.knee.x_mean;
    hip_trajectory = Gaitcycle.(subject).(trial).kinematics.jointangles.right.hip.x_mean;

    x = [hip_trajectory, knee_trajectory, ankle_trajectory];
    xstart = x(1,:);
    xstart(2) = -xstart(2);
    L0 = (fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cosd(xstart(2))+la.*ls.*cosd(xstart(3))+la.*lt.*cosd(xstart(2)+xstart(3))+(-2).*fc.*ls.*sind(xstart(3))+(-2).*fc.*lt.*sind(xstart(2)+xstart(3))).^(1/2);
    % = L0relax(j);%
    
    ToeP = [lf.*cosd(xstart(1)+xstart(2)+xstart(3))+lt.*sind(xstart(1))+ls.*sind(xstart(1)+xstart(2))+la.*sind(xstart(1)+xstart(2)+xstart(3)),(-1).*lt.*cosd(xstart(1))+(-1).*ls.*cosd(xstart(1)+xstart(2))+(-1).*la.*cosd(xstart(1)+xstart(2)+xstart(3))+lf.*sind(xstart(1)+xstart(2)+xstart(3)),0,0,0,sind(xstart(1)+xstart(2)+xstart(3))]; 
    HipP = -ToeP;
    Eref = 1/2*(1)^2*m_0 + HipP(2)*m_0*9.81;
    
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
    
    for i=1:length(ankle_trajectory)
        [Knee_torque_command(i), Ankle_torque_command(i), Knee_torque_KPBC(i), Ankle_torque_KPBC(i), deltaL(i), Ehip(i), Espring(i)]... 
        = SLIP_KPBC_forsim(hip_trajectory(i),knee_trajectory(i), ankle_trajectory(i),hip_vel(i), knee_vel(i), ankle_vel(i), lf, la, ls, lt, k, d, L0, m_0,Eref);
    end

    scale = Gaitcycle.(subject).subjectdetails{4,2}/1000; %turns N*mm/kg to N*m

    stancecutoff_value = 0.6;
    stancecutoff_index = stancecutoff_value*length(ankle_trajectory);
    
    figure('Name',subject,'NumberTitle','off')
    subplot(2,2,1)
    plot(linspace(0,stancecutoff_value,stancecutoff_index),-Ankle_torque_command(1:stancecutoff_index))
    hold on
    plot(linspace(0,stancecutoff_value,stancecutoff_index),Gaitcycle.(subject).(trial).kinetics.jointmoment.right.ankle.x_mean(1:stancecutoff_index)*scale)
    title('Ankle Torque Comparison')
    xlabel('Normalized Time')
    ylabel('Torque(N m)')
    legend('SLIP','Human')

    subplot(2,2,2)
    plot(linspace(0,stancecutoff_value,stancecutoff_index),Knee_torque_command(1:stancecutoff_index))
    hold on
    plot(linspace(0,stancecutoff_value,stancecutoff_index),Gaitcycle.(subject).(trial).kinetics.jointmoment.right.knee.x_mean(1:stancecutoff_index)*scale)
    title('Knee Torque Comparison')
    xlabel('Normalized Time')
    ylabel('Torque(N m)')
    legend('SLIP','Human')
    
    subplot(2,2,3)
    plot(linspace(0,stancecutoff_value,stancecutoff_index), ankle_trajectory(1:stancecutoff_index))
    title('Ankle Trajectory')
    xlabel('Normalized Time')
    ylabel('deg')
    
    subplot(2,2,4)
    plot(linspace(0,stancecutoff_value,stancecutoff_index), knee_trajectory(1:stancecutoff_index))
    title('Knee Trajectory')
    xlabel('Normalized Time')
    ylabel('deg')
    
%     subplot(2,2,3)
%     plot(linspace(0,stancecutoff_value,stancecutoff_index), hip_vel(1:stancecutoff_index), linspace(0,stancecutoff_value,stancecutoff_index), knee_vel(1:stancecutoff_index), linspace(0,stancecutoff_value,stancecutoff_index), ankle_vel(1:stancecutoff_index))
%     title('Velocities')
%     xlabel('Normalized Time')
%     ylabel('m/s')
%     legend('Hip','Knee','Ankle')
%     
%     subplot(2,2,4)
%     Ehip = Ehip - Ehip(1);
%     plot(linspace(0,stancecutoff_value,stancecutoff_index), Ehip(1:stancecutoff_index),linspace(0,stancecutoff_value,stancecutoff_index), Espring(1:stancecutoff_index),linspace(0,stancecutoff_value,stancecutoff_index), Ehip(1:stancecutoff_index)+ Espring(1:stancecutoff_index))
%     title('Energy')
%     xlabel('Normalized Time')
%     ylabel('J')
%     legend('E_{hip}','E_{spring}','E_{total}')
%     
%     figure('Name','Spring Length','NumberTitle','off')
%     plot(linspace(0,stancecutoff_value,stancecutoff_index),deltaL(1:stancecutoff_index))
       
end
set(0,'DefaultFigureWindowStyle','normal')

%animate(@drawProsthesis, x,[1:length(ankle_trajectory)],links)
