%istance_start and iswing_start are vectors of frames indices into a time
%vector. it will highlight stance and swing.
temp = axis;
highlight_array = {};
for i = 1:num_steps-1
    %Stance
    x = time(istance_start(i));
    w = time(iswing_start(i))-time(istance_start(i));
    y = -10000;
    h = 20000;
    pos = [x,y,w,h];
    r1 = rectangle('Position',pos,'FaceColor',[0, 1, 1, 0.3]);
    uistack(r1,'bottom')
    %alpha(r1,0.5)
    highlight_array(end+1) = {r1};
    
%     %Pushoff
%     x = time(ipush(i));
%     w = time(iswing_start(i))-time(ipush(i));
%     y = -10000;
%     h = 20000;
%     pos = [x,y,w,h];
%     r3 = rectangle('Position',pos,'FaceColor',[1, 0, 0, 0.3]);
%     uistack(r3,'bottom')
%     %alpha(r2,0.5)
%     highlight_array(end+1) = {r3};
    
    %Swing
    x = time(iswing_start(i));
    w = time(istance_start(i+1))-time(iswing_start(i));
    y = -10000;
    h = 20000;
    pos = [x,y,w,h];
    r2 = rectangle('Position',pos,'FaceColor',[1, 1, 0, 0.3]);
    uistack(r2,'bottom')
    %alpha(r2,0.5)
    highlight_array(end+1) = {r2};
end
axis(temp)
disp("Blue is stance, yellow is swing")