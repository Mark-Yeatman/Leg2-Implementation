function posAnkle = Ankle_pos_func(in1,in2)
%ANKLE_POS_FUNC
%    POSANKLE = ANKLE_POS_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    27-Jan-2020 10:54:28

%Prosthesis file. Needs state and parameters as inputs
ls = in2(:,4);
lt = in2(:,3);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
posAnkle = [x1-ls.*sin(x3);lt+x2+ls.*cos(x3);0.0;0.0;0.0;sin(x3+x4)];
