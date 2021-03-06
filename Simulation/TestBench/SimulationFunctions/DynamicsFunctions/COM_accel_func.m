function ddCoM = COM_accel_func(in1,in2,in3)
%COM_ACCEL_FUNC
%    DDCOM = COM_ACCEL_FUNC(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    31-Jan-2020 14:49:33

Mf = in3(:,3);
Ms = in3(:,2);
a1 = in2(1,:);
a2 = in2(2,:);
a3 = in2(3,:);
a4 = in2(4,:);
a5 = in2(5,:);
la = in3(:,6);
lf = in3(:,7);
ls = in3(:,5);
lt = in3(:,4);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x8 = in1(8,:);
x9 = in1(9,:);
x10 = in1(10,:);
t2 = cos(x3);
t3 = sin(x3);
t4 = a3+a4;
t5 = x3+x4;
t6 = x8+x9;
t7 = Mf.*2.0;
t8 = Mf.*4.0;
t9 = Ms.*4.0;
t10 = x8.^2;
t11 = Mf+Ms;
t12 = a5+t4;
t13 = cos(t5);
t14 = sin(t5);
t15 = t5+x5;
t16 = t6+x10;
t19 = Ms+t7;
t20 = t6.^2;
t22 = 1.0./t11;
t23 = t8+t9;
t17 = cos(t15);
t18 = sin(t15);
t21 = t16.^2;
ddCoM = [t22.*(-a1.*t23-Mf.*la.*t12.*t17.*2.0+Mf.*lf.*t12.*t18+Mf.*lf.*t17.*t21-a3.*lt.*t2.*t11.*4.0+la.*t7.*t18.*t21-ls.*t4.*t13.*t19.*2.0+ls.*t14.*t19.*t20.*2.0+lt.*t3.*t10.*t11.*4.0).*(-1.0./4.0),(t22.*(a2.*t23+Mf.*lf.*t12.*t17-Mf.*lf.*t18.*t21+a3.*lt.*t3.*t11.*4.0+la.*t7.*t12.*t18+la.*t7.*t17.*t21+ls.*t4.*t14.*t19.*2.0+ls.*t13.*t19.*t20.*2.0+lt.*t2.*t10.*t11.*4.0))./4.0,0.0];
