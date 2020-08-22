function Pe = PE_func(in1,in2)
%PE_FUNC
%    PE = PE_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    27-Jan-2020 14:19:42

%Prosthesis file. Needs state and parameters as inputs
Mf = in2(:,2);
Ms = in2(:,1);
la = in2(:,5);
lf = in2(:,6);
ls = in2(:,4);
lt = in2(:,3);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x3+x4;
Pe = x2.*(Mf.*4.0+Ms.*4.0).*(9.81e+2./4.0e+2)-lt.*(Mf+Ms).*(9.81e+2./1.0e+2)+Mf.*lf.*sin(t2).*(9.81e+2./4.0e+2)-ls.*cos(x3).*(Mf.*2.0+Ms).*(9.81e+2./2.0e+2)-Mf.*la.*cos(t2).*(9.81e+2./2.0e+2);
end