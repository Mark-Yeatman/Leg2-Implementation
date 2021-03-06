function lJacob = Spring_Jacobian_func(in1,in2)
%SPRING_JACOBIAN_FUNC
%    LJACOB = SPRING_JACOBIAN_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    31-Jan-2020 14:49:45

%Prosthesis file. Needs state and parameters as inputs
la = in2(:,6);
lf = in2(:,7);
ls = in2(:,5);
lt = in2(:,4);
px = in2(:,8);
py = in2(:,9);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
t2 = cos(x3);
t3 = sin(x3);
t4 = la.*py;
t5 = lf.*px;
t6 = x3+x4;
t7 = x4+x5;
t15 = -px;
t8 = lt.*t2;
t9 = cos(t6);
t10 = cos(t7);
t11 = lt.*t3;
t12 = sin(t6);
t13 = sin(t7);
t14 = t6+x5;
t18 = -t5;
t16 = cos(t14);
t17 = sin(t14);
t19 = ls.*t9;
t20 = ls.*t12;
t23 = lf.*lt.*t10;
t27 = la.*lt.*t13;
t32 = t4+t18;
t21 = la.*t17;
t22 = lf.*t17;
t24 = px.*t19;
t25 = la.*t16;
t26 = lf.*t16;
t28 = py.*t20;
t33 = t17.*t32;
t29 = px.*t25;
t30 = py.*t26;
t31 = -t22;
t35 = t11+t15+t20+t21+t26;
t34 = py+t8+t19+t25+t31;
t37 = t35.^2;
t36 = t34.^2;
t38 = t36+t37;
t39 = 1.0./sqrt(t38);
lJacob = [0.0;0.0;-t39.*(t24+t28+t29+t30+t33+px.*t8+py.*t11);-t39.*(t23+t24+t27+t28+t33+t16.*(la.*px+lf.*py)+ls.*lt.*sin(x4));-t39.*(t23+t27+t29+t30+t33+lf.*ls.*cos(x5)+la.*ls.*sin(x5))];
end