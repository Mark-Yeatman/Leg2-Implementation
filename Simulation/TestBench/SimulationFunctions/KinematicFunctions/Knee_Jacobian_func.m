function js0 = Knee_Jacobian_func(in1,in2)
%KNEE_JACOBIAN_FUNC
%    JS0 = KNEE_JACOBIAN_FUNC(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    31-Jan-2020 14:49:44

%Prosthesis file. Needs state and parameters as inputs
lt = in2(:,4);
x3 = in1(3,:);
t2 = cos(x3);
js0 = reshape([1.0,1.0,0.0,0.0,0.0,0.0,0.0,lt.*sin(x3),0.0,0.0,0.0,0.0,lt.*t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,6]);
