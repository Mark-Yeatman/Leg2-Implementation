figure('Name','Step Info','NumberTitle','off');

axSI1 = subplot(3,1,1);
plot(time,foot_contact(I))
xlabel("Time (s)")
title("Foot Contact")

axSI2 = subplot(3,1,2);
plot(time,PushOff(I))
xlabel("Time (s)")
title("PushOff");

axSI3 = subplot(3,1,3);
plot(time, hip_pos(I));
xlabel("Time (s)")
ylabel("Angle")
title("Hip Position")

linkaxes([axSI1, axSI2, axSI3],'x')