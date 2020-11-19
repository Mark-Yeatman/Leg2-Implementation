function posToe = Toe_pos_func(in1,in2)
%TOE_POS_FUNC
%    POSTOE = TOE_POS_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    27-Jan-2020 10:54:28

%Prosthesis file. Needs state and parameters as inputs
la = in2(:,5);
lf = in2(:,6);
ls = in2(:,4);
lt = in2(:,3);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x3+x4;
t3 = cos(t2);
t4 = sin(t2);
posToe = [x1-la.*t4+lf.*t3-ls.*sin(x3);lt+x2+la.*t3+lf.*t4+ls.*cos(x3);0.0;0.0;0.0;t4];