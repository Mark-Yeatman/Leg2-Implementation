PushOff = zeros(size(time_in));
for i = 1:length(time_in)
    PushOff(i) = getPush(hip_pos(i),hip_thresh_st(i),stance(i));
end

function PushOffOut = getPush(hip_pos, hip_thresh_st, Stance)
    %GETPUSHOFF Summary of this function goes here
    %   Detailed explanation goes here
    persistent PushOff
    if isempty(PushOff)
       PushOff = false; 
    end
    PushOff = and(or((hip_pos < hip_thresh_st),PushOff),Stance);
    PushOffOut = PushOff;
end

