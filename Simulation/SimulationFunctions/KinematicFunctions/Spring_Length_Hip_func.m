function LFlatHip = Spring_Length_Hip_func(in1)
%SPRING_LENGTH_HIP_FUNC
%    LFLATHIP = SPRING_LENGTH_HIP_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    14-Oct-2019 12:19:56

%Prosthesis file. Needs state and parameters as inputs
global flowdata
in2 =  flowdata.Parameters.Biped.asvector;
la = in2(:,10);
lf = in2(:,11);
ls = in2(:,9);
lt = in2(:,8);
y2 = in1(2,:);
y3 = in1(3,:);
t2 = y2+y3;
LFlatHip = sqrt(la.^2+lf.^2+ls.^2+lt.^2+la.*lt.*cos(t2).*2.0+la.*ls.*cos(y3).*2.0+ls.*lt.*cos(y2).*2.0-lf.*lt.*sin(t2).*2.0-lf.*ls.*sin(y3).*2.0);
