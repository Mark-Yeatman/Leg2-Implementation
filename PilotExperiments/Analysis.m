%% Time window, ts_i and ts_e should be set manually
[time,I] = sort(time_in-min(time_in));

if exist('ts_i','var') && exist('te_i','var')
    time = time(ts_i:te_i);
    I = I(ts_i:te_i);
else
    disp("No ts_i,te_i variables, starting from t = 0")
end

%% Stats

%Find heel strikes \ steps
dfc = gradient(foot_contact(I), time);                
[mag,istance_start] = findpeaks(dfc, 'MinPeakDistance',0.1, 'MinPeakHeight',0.5);
[~,iswing_start] = findpeaks(-dfc, 'MinPeakDistance',0.1, 'MinPeakHeight',0.5);
if iswing_start(1)<istance_start(1)
    iswing_start = iswing_start(2:end); %cut down swing start indices to start with stance
end

percent_gait = zeros(size(time));
for i = 1:length(time)
    t_start = time(istance_start( find(time(istance_start)<=time(i),1,'last') ));
    t_end = time(istance_start( find(time(istance_start)> time(i),1,'first') ));
    if isempty(t_start) || isempty(t_end)
        percent_gait(i) = nan;
    else
        percent_gait(i) = (time(i)-t_start)/(t_end-t_start);
    end
end

num_steps = length(istance_start);
sprintf("Number of steps: %i",num_steps) 

%Find outlier steps
steptimes = diff(time(istance_start)); %doesn't count last step
outlier_steps = isoutlier(steptimes,'grubbs');
disp("Outlier steps by step time: ")
disp(find(outlier_steps==true))

%getJointVelocities;
%PushOff = (hip_pos(I) > hip_thresh_st(I)) & (hip_vel(I)>hip_vel_thresh_st(I));
%dfc = gradient(double(PushOff), time);                
%[mag,ipush] = findpeaks(dfc, 'MinPeakDistance',0.1, 'MinPeakHeight',0.5);