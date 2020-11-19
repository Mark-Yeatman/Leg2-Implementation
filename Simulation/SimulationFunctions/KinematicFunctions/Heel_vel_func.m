function velHeel = Heel_vel_func(in1,in2)
%HEEL_VEL_FUNC
%    VELHEEL = HEEL_VEL_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    27-Jan-2020 10:54:29

%Prosthesis file. Needs state and parameters as inputs
la = in2(:,5);
ls = in2(:,4);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
x7 = in1(7,:);
x8 = in1(8,:);
t2 = x3+x4;
t3 = x7+x8;
t4 = cos(t2);
velHeel = [x5-ls.*x7.*cos(x3)-la.*t3.*t4;x6-la.*t3.*sin(t2)-ls.*x7.*sin(x3);0.0;0.0;0.0;t3.*t4];
