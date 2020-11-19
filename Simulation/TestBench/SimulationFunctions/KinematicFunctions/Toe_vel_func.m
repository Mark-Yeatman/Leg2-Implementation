function velToe = Toe_vel_func(in1,in2)
%TOE_VEL_FUNC
%    VELTOE = TOE_VEL_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    31-Jan-2020 14:49:43

%Prosthesis file. Needs state and parameters as inputs
la = in2(:,6);
lf = in2(:,7);
ls = in2(:,5);
lt = in2(:,4);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
x7 = in1(7,:);
x8 = in1(8,:);
x9 = in1(9,:);
x10 = in1(10,:);
t2 = x3+x4;
t3 = x8+x9;
t4 = t2+x5;
t5 = t3+x10;
t6 = cos(t4);
t7 = sin(t4);
velToe = [x6+ls.*t3.*cos(t2)+lt.*x8.*cos(x3)+la.*t5.*t6-lf.*t5.*t7;x7+ls.*t3.*sin(t2)+lt.*x8.*sin(x3)+la.*t5.*t7+lf.*t5.*t6;0.0;0.0;0.0;t5.*t6];
