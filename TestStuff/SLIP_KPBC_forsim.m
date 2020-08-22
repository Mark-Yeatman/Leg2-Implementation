function [Knee_torque_slip, Ankle_torque_slip, Knee_torque_KPBC, Ankle_torque_KPBC, deltaL, Ehip, Espring]...
= SLIP_KPBC_forsim(hip_pos, knee_pos, ankle_pos, hip_vel, knee_vel, ankle_vel, lf, la, ls, lt, k, d, L0, M, Eref)
   
    %Assign joint positions
    knee_pos = - knee_pos; %reverse knee sign convention of biomechanics versus biped modeling
    
    x = zeros(6,1);
    x(1) =  hip_pos;
    x(2) =  knee_pos;         
    x(3) =  ankle_pos;
    x(4) =  hip_vel;
    x(5) =  -knee_vel; %because its not calculated using a finite difference method, its just given to us right now, change for implementation
    x(6) =  ankle_vel;
     
    %% Embedded SLIP Controller %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fc = lf;
    
    J = [0;(-1).*lt.*(2.*fc.*cosd(x(2)+x(3))+2.*ls.*sind(x(2))+la.*sind(x(2)+x(3))).*(4.*fc.^2+la.^2+4.*ls.^2+4.*lt.^2+8.*ls.*lt.*cosd(x(2))+4.*la.*ls.*cosd(x(3))+4.*la.*lt.*cosd(x(2)+x(3))+(-8).*fc.*ls.*sind(x(3))+(-8).*fc.*lt.*sind(x(2)+x(3))).^(-1/2);(-1).*(4.*fc.^2+la.^2+4.*ls.^2+4.*lt.^2+8.*ls.*lt.*cosd(x(2))+4.*la.*ls.*cosd(x(3))+4.*la.*lt.*cosd(x(2)+x(3))+(-8).*fc.*ls.*sind(x(3))+(-8).*fc.*lt.*sind(x(2)+x(3))).^(-1/2).*(2.*fc.*ls.*cosd(x(3))+2.*fc.*lt.*cosd(x(2)+x(3))+la.*ls.*sind(x(3))+la.*lt.*sind(x(2)+x(3)))];
    L  = (fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cosd(x(2))+la.*ls.*cosd(x(3))+la.*lt.*cosd(x(2)+x(3))+(-2).*fc.*ls.*sind(x(3))+(-2).*fc.*lt.*sind(x(2)+x(3))).^(1/2);
    Ldot = (1/2).*(fc.^2+(1/4).*la.^2+ls.^2+lt.^2+2.*ls.*lt.*cos(x(2))+la.*ls.*cos(x(3))+la.*lt.*cos(x(2)+x(3))+(-2).*fc.*ls.*sin(x(3))+(-2).*fc.*lt.*sin(x(2)+x(3))).^(-1/2).*((-2).*ls.*lt.*sin(x(2)).*x(5)+(-2).*fc.*ls.*cos(x(3)).*x(6)+(-1).*la.*ls.*sin(x(3)).*x(6)+(-2).*fc.*lt.*cos(x(2)+x(3)).*(x(5)+x(6))+(-1).*la.*lt.*sin(x(2)+x(3)).*(x(5)+x(6)));
    deltaL = L - L0;
    u = -k.*(deltaL).*J - d.*(Ldot).*J;
    Espring =  1/2*k*(deltaL^2);
    Knee_torque_slip  = u(2);
    Ankle_torque_slip = u(3);

    %% Energy Tracking Controller %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = x(4:6);
    ToeP = [lf.*cosd(x(1)+x(2)+x(3))+lt.*sind(x(1))+ls.*sind(x(1)+x(2))+la.*sind(x(1)+x(2)+x(3)),(-1).*lt.*cosd(x(1))+(-1).*ls.*cosd(x(1)+x(2))+(-1).*la.*cosd(x(1)+x(2)+x(3))+lf.*sind(x(1)+x(2)+x(3)),0,0,0,sind(x(1)+x(2)+x(3))]; 
    HipP = -ToeP; %because of the switch from HipFrame to ToeFrame
    ToeV = [-(v(1)+v(2)+v(3))*lf.*sind(x(1)+x(2)+x(3)) + v(1)*lt.*cosd(x(1)) + (v(1)+v(2))*ls.*cosd(x(1)+x(2)) + (v(1)+v(2)+v(3))*la.*cosd(x(1)+x(2)+x(3)),...
       v(1).*lt.*sind(x(1)) + (v(1)+v(2)).*ls.*sind(x(1)+x(2)) + (v(1)+v(2)+v(3)).*la.*sind(x(1)+x(2)+x(3)) + (v(1)+v(2)+v(3))*lf.*cosd(x(1)+x(2)+x(3)),...
       0,...
       0,...
       0,...
       (v(1)+v(2)+v(3))*cosd(x(1)+x(2)+x(3))]; 
    HipV = -ToeV; %because of the switch from HipFrame to ToeFrame

    Ehip = norm(HipV(1:2),2)*M*1/2 + HipP(2)*M*9.81;%HipP(2);%
    y = -Ldot*(Ehip-Eref);
    u = y*J;
    Knee_torque_KPBC  = u(2);
    Ankle_torque_KPBC = u(3);
       

end
