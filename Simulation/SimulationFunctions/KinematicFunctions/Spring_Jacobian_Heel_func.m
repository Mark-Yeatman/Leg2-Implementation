function LJacobFlatHeel = Spring_Jacobian_Heel_func(in1)
%SPRING_JACOBIAN_HEEL_FUNC
%    LJACOBFLATHEEL = SPRING_JACOBIAN_HEEL_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    14-Oct-2019 12:19:58

%Prosthesis file. Needs state and parameters as inputs
global flowdata
in2 =  flowdata.Parameters.Biped.asvector;
la = in2(:,10);
lf = in2(:,11);
ls = in2(:,9);
lt = in2(:,8);
x4 = in1(4,:);
x5 = in1(5,:);
t2 = cos(x4);
t3 = cos(x5);
t4 = sin(x4);
t5 = x4+x5;
t6 = la.^2;
t7 = lf.^2;
t8 = ls.^2;
t9 = lt.^2;
t10 = cos(t5);
t11 = sin(t5);
t12 = la.*ls.*t2.*2.0;
t13 = ls.*lt.*t3.*2.0;
t14 = lf.*ls.*t4.*2.0;
t15 = la.*lt.*t10.*2.0;
t16 = lf.*lt.*t11.*2.0;
t17 = t6+t7+t8+t9+t12+t13+t14+t15+t16;
t18 = 1.0./sqrt(t17);
LJacobFlatHeel = [0.0;0.0;0.0;t18.*(-la.*(ls.*t4+lt.*t11)+lf.*ls.*t2+lf.*lt.*t10);-lt.*t18.*(la.*t11-lf.*t10+ls.*sin(x5));0.0;0.0;0.0];