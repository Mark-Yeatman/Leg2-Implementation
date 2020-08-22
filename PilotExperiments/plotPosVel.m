figure('Name','Pos and Vel','NumberTitle','off');

tempKnee = knee_des_out(I);
tempKnee(foot_contact(I) == 1) = nan;
tempAnkle = ankle_des_out(I);
tempAnkle(foot_contact(I) == 1) = nan;
ax5 = subplot(2,3,1);
plot(time, knee_joint_pos(I), time, ankle_joint_pos(I), time, tempKnee, time, tempAnkle);
xlabel("Time (s)")
ylabel("Angle (deg)")
legend("Knee Pos","Ankle Pos","Desired Knee", "Desired Ankle")

ax6 = subplot(2,3,2);
plot(time, knee_joint_vel(I),time,ankle_joint_vel(I));
legend("Knee Vel","Ankle Vel")
xlabel("Time (s)")
ylabel("Angular Velocity (deg/s)")

ax36 = subplot(2,3,3);
plot(time, hip_pos(I));
xlabel("Time (s)")
ylabel("Angle")
title("Hip Position")

ax37 = subplot(2,3,4);
plot(time, phase_out(I));
xlabel("Time (s)")
ylabel("Phase ")
title("Phase Variable")

ax38 = subplot(2,3,5);
plot(time, delta_L(I));
xlabel("Time (s)")
ylabel("Dist (m)")
title("Virtual Spring Displacement")

Ldot = 0*delta_L;
for i=1:length(Ldot)
    x = zeros(10,1); %x,y,-knee ankle, ankle angle
    x(4) =  -knee_joint_pos(i); %reverse knee sign convention of biomechanics versus biped modeling     
    x(5) =  ankle_joint_pos(i);

    x(9) =  -knee_joint_vel(i); %reverse knee sign convention of biomechanics versus biped modeling
    x(10) =  ankle_joint_vel(i);
    Mt = 0;
    Ms = single(6.991429/2.205); % converted to kg
    Mf = single(0.55 + 1.643130/2.205); %carbon fiber foot + mechanism, kg
    %lt = single(0.3733); %meters
    ls = single(0.3733); %meters
    la = single(0.0628); %meters
    lf = single(0.15); %meters
    px = single(0); %meters
    py = single(0); %meters
    params = [Mt, Ms, Mf, lt1(i), ls, la, lf, px, py]; %From ordering in makeMatlabFunctionsProthesisTestBench
    Ldot(i) = Spring_vel_func(x,params);
end

ax39 = subplot(2,3,6);
plot(time, Ldot(I));
xlabel("Time (s)")
ylabel("Vel (m/s)")
title("Virtual Spring Velocity")

linkaxes([ax5,ax6, ax36, ax37,ax38, ax39],'x');