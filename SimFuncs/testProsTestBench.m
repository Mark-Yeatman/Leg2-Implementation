%clear all
path(pathdef)
addpath('Experiments\ModelTests\')
addpath('Analysis\')
addpath('UtilityFunctions\')

%Get all the simulations functions
addpath(genpath('Models\Prosthesis\TestBench'))

global flowdata
flowdata = flowData;

flowdata.E_func = @E_func;

%Ode equation handle and tolerenaces
flowdata.eqnhandle = @dynamics;
flowdata.odeoptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%Flags
flowdata.Flags.silent = false;
flowdata.Flags.ignore = true;

%Simulation parameters
flowdata.Parameters.Environment.slope = 0;   %ground slope
dim = 8;
flowdata.Parameters.dim = dim;                 %state variable dimension

%Dynamics Parameters, from Leg 2 Overview.docx
Ms = 6.991429/2.205; % converted to kg
Mf = 0.55 + 1.643130/2.205; %carbon fiber foot + mechanism, kg
lt = 0.0959; %meters
ls = 0.3733; %meters
la = 0.0628; %meters
lf = 0.15; %meters
px = 0; %meters
py = 0; %meters
flowdata.Parameters.Dynamics.asvector = [Ms, Mf, lt, ls, la, lf, px, py]; %From ordering in makeMatlabFunctionsProthesisTestBench
%px, py is hip side spring mount point, other is the toe

%Discrete Mappings and Constraints
flowdata.setPhases({})
flowdata.setConfigs({'Mounted'})
impactlist =  {};
n_phaselist = {};
n_configlist = {'Mounted'};
flowdata.setImpacts(impactlist,n_phaselist,n_configlist);
flowdata.End_Step.event_name = '';

%Set initial phase and contact conditions
flowdata.State.c_phase = {};
flowdata.State.c_configs = {'Mounted'};
flowdata.State.Eref = 0;
flowdata.odeoptions.Events = [];
flowdata.tspan = 5; %5 seconds of in simulation time

%Load initial condition 
xi = zeros(1,dim);
x_E_min = [0, 0,-0.000000052891209, -0.873702184821549 ,0, 0 ,0,0];
perturb = [0, 0,-0.1, -0.1 ,0, 0 ,0,0];
E_min = -15.371619960146401;
%Control Functions
%flowdata.Controls.Internal = {@SLIP, @KPBC};
%flowdata.Controls.Internal = {@PD};
flowdata.Controls.Internal = {@SLIP,@KPBC};

%Control Parameters
flowdata.Parameters.PD.KD = [0.25,0.5];
flowdata.Parameters.PD.KP = 0;
flowdata.Parameters.PD.setpoint = [0;0];

flowdata.Parameters.SLIP.L0 = 0.839381199857126;
flowdata.Parameters.SLIP.k = 10; 
flowdata.Parameters.SLIP.d = 0;

flowdata.Parameters.KPBC.k = 1;
flowdata.Parameters.KPBC.omega = diag([0,0,1,0.01]);
flowdata.Parameters.KPBC.sat = 100;
flowdata.Parameters.KPBC.Eref = 0.1+E_min;

%flowdata.E_func = @Energy_Sub;
flowdata.Flags.do_validation = false;

flowdata.tspan = 0.25;
[fstate, xout, tout, out_extra] = walk(x_E_min+perturb,1);

videopath = 'Experiments\Videos\';
animate(@drawProsthesisTestBench,xout,tout,[],1,strcat(videopath,'testProsTestBench_OpenLoop'))
%animate(@drawProsthesisTestBench,xout,tout,[],1)