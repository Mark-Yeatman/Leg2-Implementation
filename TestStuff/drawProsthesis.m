function drawProsthesis(x, links)
    %DRAWPROSTHESIS Draw picture of hip/thigh wearing a trans femural
    %prosthesis.
    %   Detailed explanation goes here
    
    lt = links(1);
    ls = links(2);
    la = links(3);
    lf = links(4);
    
    x(2) = -x(2);
    
    Knee = [lt.*sind(x(1)),(-1).*lt.*cosd(x(1)),0,0,0,sind(x(1))];
    Ankle = [lt.*sind(x(1))+ls.*sind(x(1)+x(2)),(-1).*lt.*cosd(x(1))+(-1).*ls.*cosd(x(1)+x(2)),0,0,0,sind(x(1)+x(2))];
    Heel = [lt.*sind(x(1))+ls.*sind(x(1)+x(2))+la.*sind(x(1)+x(2)+x(3)),(-1).*lt.*cosd(x(1))+(-1).*ls.*cosd(x(1)+x(2))+(-1).*la.*cosd(x(1)+x(2)+x(3)),0,0,0,sind(x(1)+x(2)+x(3))];
    Toe = [lf.*cosd(x(1)+x(2)+x(3))+lt.*sind(x(1))+ls.*sind(x(1)+x(2))+la.*sind(x(1)+x(2)+x(3)),(-1).*lt.*cosd(x(1))+(-1).*ls.*cosd(x(1)+x(2))+(-1).*la.*cosd(x(1)+x(2)+x(3))+lf.*sind(x(1)+x(2)+x(3)),0,0,0,sind(x(1)+x(2)+x(3))];
    limb_x = [0, Knee(1)];
    limb_y = [0, Knee(2)];
    pros_x = [Knee(1), Ankle(1), Heel(1), Toe(1)];
    pros_y = [Knee(2), Ankle(2), Heel(2), Toe(2)];
    
    hold off %clears previous content with next plot command
    plot(limb_x, limb_y, 'Color','g','LineWidth',3)
    hold on
    plot(pros_x, pros_y,'k', 'LineWidth',3)    
    
    %Draw joint points
    plot([limb_x, pros_x], [limb_y, pros_y], 'o', 'LineWidth',1,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
    
    axis equal
    axis([-1,1,-1,1])
    grid on
end

