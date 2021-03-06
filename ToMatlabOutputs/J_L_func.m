function JacobL = J_L_func(in1,in2)
%J_L_FUNC
%    JACOBL = J_L_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    25-Aug-2020 13:50:39

COPfx = in2(:,12);
la = in2(:,8);
ls = in2(:,6);
lt = in2(:,5);
x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1+x2;
t3 = -lt;
t5 = (x1.*pi)./1.8e+2;
t6 = (x2.*pi)./1.8e+2;
t4 = ls+t3;
t7 = cos(t5);
t8 = sin(t5);
t9 = (t2.*pi)./1.8e+2;
t10 = cos(t9);
t11 = sin(t9);
t17 = t4.*t7;
t18 = t4.*t8;
t12 = COPfx.*t10;
t13 = COPfx.*t11;
t14 = la.*t10;
t15 = la.*t11;
t16 = -t13;
t19 = t12+t15+t18;
t20 = t19.^2;
t21 = lt+t14+t16+t17;
t22 = t21.^2;
t23 = t20+t22;
t24 = 1.0./sqrt(t23);
JacobL = [t3.*t19.*t24,-t24.*(la.*(lt.*t11+t4.*sin(t6))+lt.*t12+COPfx.*t4.*cos(t6))];
