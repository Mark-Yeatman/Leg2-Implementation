winter_data_path = '\\engstor\locolab\Tutorials\Human Gait\winter_dataset_std.mat';
load(winter_data_path)

figure('Name','Winter Data Comparison','NumberTitle','off');

%In Winter's 'Biomechanics and Motor Control of Human Movement', pg 194,
%the sign convention for torque is backwards compared to the convention for
%position
axWC1 = subplot(2,1,1);
kneeW_data = level_walking_std.normal_cadence.knee.position;
kneeW_mean = interp1(kneeW_data.stride_percent/100,kneeW_data.avg,percent_gait);
kneeW_std = interp1(kneeW_data.stride_percent/100,kneeW_data.std,percent_gait);
kneeW_std_low = kneeW_mean - kneeW_std;
kneeW_std_high = kneeW_mean + kneeW_std;

plot(time,knee_joint_pos,time,kneeW_mean)
xlabel('Time (s)')
ylabel('Pos (degrees)')
title('Knee')
legend('Pros','Winters')
%hold on 
%plot(time,kneeW_mean,'b',time,kneeW_mean,'b')

axWC2 = subplot(2,1,2);
ankleW_data = level_walking_std.normal_cadence.ankle.position;
ankleW_mean = interp1(ankleW_data.stride_percent/100,ankleW_data.avg,percent_gait);
ankleW_std = interp1(ankleW_data.stride_percent/100,ankleW_data.std,percent_gait);
ankleW_std_low = ankleW_mean - ankleW_std;
ankleW_std_high = ankleW_mean + ankleW_std;

plot(time,ankle_joint_pos,time,ankleW_mean)
xlabel('Time (s)')
ylabel('Pos (degrees)')
title('Ankle')
legend('Pros','Winters')
%hold on 
%plot(time,ankleW_mean,'b',time,ankleW_mean,'b')