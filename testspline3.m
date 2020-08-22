function [y,dy]=testspline3(x,y0,dy0,y1,dy1)
%x is spline input parameter between 0 and 1
% y0 and y1 are the position start and finish
% dy0 and dy1 are the velocity start and finish
    a0=y0;
    a1=dy0;
    a2=3*y1-3*y0-2*dy0-dy1;
    a3=-2*y1+2*y0+dy0+dy1;
    if x<0
        y=y0;
        dy=dy0;
    elseif x>1
        y=y1;
        dy=dy1;
    else
        y=a0+a1*x+a2*x.^2+a3*x.^3;
        dy=a1+2*a2*x+3*a3*x.^2;
    end
end