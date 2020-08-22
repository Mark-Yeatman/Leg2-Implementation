%Notes:
% Knee axis treated as origin.
% Orientation of knee and ankle rotation is biomechanical

%variables:
% x - 4x1 array [knee pos, ankle pos, knee vel, ankle vel]
% lt -thigh length 
% ls -shank length
% lfx, lfy - foot CoM x,y
% lc - load cell y dist from ankle axis
% la - ankle axis y dist to bottom of foot
% COPfx - COP distance along foot from ankle axis projection, to be
%   computed from load cell
% Ms - shank mass
% Mf - foot mass

M=[(1/4).*lt.^2.*Ms+Mf.*(lfy+(-1).*ls+lt+ls.*cos(x(2))).^2+Mf.*(lfx+(-1).*ls.*sin(x(2))).^2,(lfy+(-1).*ls+lt).*Mf.*(lfy+(-1).*ls+lt+ls.*cos(x(2)))+lfx.*Mf.*(lfx+(-1).*ls.*sin(x(2)));(lfy+(-1).*ls+lt).*Mf.*(lfy+(-1).*ls+lt+ls.*cos(x(2)))+lfx.*Mf.*(lfx+(-1).*ls.*sin(x(2))),lfx.^2.*Mf+(lfy+(-1).*ls+lt).^2.*Mf];
C=[ls.*Mf.*((-1).*lfx.*cos(x(2))+(-1).*(lfy+(-1).*ls+lt).*sin(x(2))).*x(4),ls.*Mf.*((-1).*lfx.*cos(x(2))+(-1).*(lfy+(-1).*ls+lt).*sin(x(2))).*(x(3)+x(4));ls.*Mf.*(lfx.*cos(x(2))+(lfy+(-1).*ls+lt).*sin(x(2))).*x(3),0];
G=[g.*(lfx.*Mf.*cos(x(1)+x(2))+(ls.*Mf+(1/2).*lt.*Ms).*sin(x(1))+(lfy+(-1).*ls+lt).*Mf.*sin(x(1)+x(2)));g.*(lfx.*Mf.*cos(x(1)+x(2))+(lfy+(-1).*ls+lt).*Mf.*sin(x(1)+x(2)))];
JC=[lc+(-1).*ls+lt+ls.*cos(x(2)),(-1).*ls.*sin(x(2)),0,0,0,1;lc+(-1).*ls+lt,0,0,0,0,1];

% COPfx calc
n = [0,1,0]';
MA = [Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment]';
F = [Load_cell_x_force, Load_cell_y_force, Load_cell_z_force]';

OA = la*n;
M0 = (cross(OA,F)+MA);

OC = cross(n,M0)/dot(F,n);
COPfx = OC(1);

L=((lt+ls.*cos(x(1))+la.*cos(x(1)+x(2))+(-1).*COPfx.*sin(x(1)+x(2))).^2+(COPfx.*cos(x(1)+x(2))+ls.*sin(x(1))+la.*sin(x(1)+x(2))).^2).^(1/2);
JL=[(-1).*lt.*(COPfx.*cos(x(1)+x(2))+ls.*sin(x(1))+la.*sin(x(1)+x(2))).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).^(-1/2);(-1).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).^(-1/2).*(COPfx.*ls.*cos(x(2))+COPfx.*lt.*cos(x(1)+x(2))+la.*ls.*sin(x(2))+la.*lt.*sin(x(1)+x(2)))];
HL=[(1/2).*lt.*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).^(-3/2).*((-2).*lt.*(COPfx.*cos(x(1)+x(2))+ls.*sin(x(1))+la.*sin(x(1)+x(2))).^2+(-2).*(ls.*cos(x(1))+la.*cos(x(1)+x(2))+(-1).*COPfx.*sin(x(1)+x(2))).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2)))),(1/2).*lt.*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).^(-3/2).*(2.*((-1).*la.*cos(x(1)+x(2))+COPfx.*sin(x(1)+x(2))).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2)))+(-2).*(COPfx.*cos(x(1)+x(2))+ls.*sin(x(1))+la.*sin(x(1)+x(2))).*(COPfx.*ls.*cos(x(2))+COPfx.*lt.*cos(x(1)+x(2))+la.*ls.*sin(x(2))+la.*lt.*sin(x(1)+x(2))));(1/2).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).^(-3/2).*(2.*lt.*((-1).*la.*cos(x(1)+x(2))+COPfx.*sin(x(1)+x(2))).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2)))+(-2).*lt.*(COPfx.*cos(x(1)+x(2))+ls.*sin(x(1))+la.*sin(x(1)+x(2))).*(COPfx.*ls.*cos(x(2))+COPfx.*lt.*cos(x(1)+x(2))+la.*ls.*sin(x(2))+la.*lt.*sin(x(1)+x(2)))),(1/2).*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).^(-3/2).*(2.*(COPfx.^2+la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(1))+2.*la.*ls.*cos(x(2))+2.*la.*lt.*cos(x(1)+x(2))+(-2).*COPfx.*ls.*sin(x(2))+(-2).*COPfx.*lt.*sin(x(1)+x(2))).*((-1).*la.*ls.*cos(x(2))+(-1).*la.*lt.*cos(x(1)+x(2))+COPfx.*ls.*sin(x(2))+COPfx.*lt.*sin(x(1)+x(2)))+(-2).*(COPfx.*ls.*cos(x(2))+COPfx.*lt.*cos(x(1)+x(2))+la.*ls.*sin(x(2))+la.*lt.*sin(x(1)+x(2))).^2)];
