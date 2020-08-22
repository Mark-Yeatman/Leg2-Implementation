winter_data_path = "Z:\Tutorials\Human Gait\winter_dataset_std.mat";
load(winter_data_path)

figure('Name','Winter Data Comparison','NumberTitle','off');

%In Winter's "Biomechanics and Motor Control of Human Movement", pg 194,
%the sign convention for torque is backwards compared to the convention for
%position
mass = 91; %kg
axWC1 = subplot(2,1,1);
kneeW_data = level_walking_std.normal_cadence.knee.torque;
kneeW_mean = -interp1(kneeW_data.stride_percent/100,kneeW_data.avg*mass,percent_gait);
kneeW_std = -interp1(kneeW_data.stride_percent/100,kneeW_data.std*mass,percent_gait);
kneeW_std_low = kneeW_mean - kneeW_std;
kneeW_std_high = kneeW_mean + kneeW_std;

plot(time,Knee_torque_command(I),time,kneeW_mean)
xlabel("Time (s)")
ylabel("Torque (N m)")
title("Knee")
legend("Pros","Winters")
%hold on 
%plot(time,kneeW_mean,"b",time,kneeW_mean,"b")

axWC2 = subplot(2,1,2);
ankleW_data = level_walking_std.normal_cadence.ankle.torque;
ankleW_mean = -interp1(ankleW_data.stride_percent/100,ankleW_data.avg*mass,percent_gait);
ankleW_std = -interp1(ankleW_data.stride_percent/100,ankleW_data.std*mass,percent_gait);
ankleW_std_low = ankleW_mean - ankleW_std;
ankleW_std_high = ankleW_mean + ankleW_std;

plot(time,Ankle_torque_command(I),time,ankleW_mean)
xlabel("Time (s)")
ylabel("Torque (N m)")
title("Ankle")
legend("Pros","Winters")
%hold on 
%plot(time,ankleW_mean,"b",time,ankleW_mean,"b")