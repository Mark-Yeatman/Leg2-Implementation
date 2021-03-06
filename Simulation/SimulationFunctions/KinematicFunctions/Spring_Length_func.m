function l1 = Spring_Length_func(in1,in2)
%SPRING_LENGTH_FUNC
%    L1 = SPRING_LENGTH_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    27-Jan-2020 10:54:28

%Prosthesis file. Needs state and parameters as inputs
la = in2(:,5);
lf = in2(:,6);
ls = in2(:,4);
lt = in2(:,3);
px = in2(:,7);
py = in2(:,8);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x3+x4;
t3 = la.*2.0;
t6 = -lt;
t4 = cos(t2);
t5 = sin(t2);
t7 = t3+t6;
l1 = sqrt((lt-py+lf.*t5-t4.*(lt-t3)+ls.*cos(x3)).^2+(px-lf.*t4-t5.*(lt-t3)+ls.*sin(x3)).^2);
