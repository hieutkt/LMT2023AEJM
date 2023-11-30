% model_baseline_PE.m
%
%   Baseline model in PE.
%
%   1. Uses stationary equilibrium solutions as beginning and end points.
%   2. Uses {P(t)} from GE solution while imposing fixed interest rates.
%
function model_baseline_PE(model)

%% Housekeeping
clc
close all

% set path to send outputs to / xstar location
mat_out = 'baseline_PE/mat_files';
mkdir(mat_out);

% fixed R (PE) or vary R (GE)?
par_in.fixedR = 1;

par_in.new_Pnow = 1/3; % fixed P (Clementi+Palazzo calibration)
par_in.Q = par_in.new_Pnow;


%% Run the steady states and transitions

%%% solve transition path
t_trans.Pf = 0;
t_trans.T = 1;
t_trans.Tmax = 40;
t_trans.fixedR = par_in.fixedR;
t_trans.Q = par_in.Q;
% path to baseline GE solution
file_in = 'mat_files_compiled/baseline/mat_files/';
switch model
    case 'baseline'
        load('mat_files_compiled/baseline/mat_files/transition_baseline.mat','P_t')
        exovars.P_t = P_t;
        
        f_transition_PE_convex('ss_autarky_baseline','ss_tradeshock_baseline',file_in,mat_out,t_trans,'_PE',exovars)
        
    case 'frictionless'
        load('mat_files_compiled/baseline/mat_files/transition_q1_baseline.mat','P_t')
        exovars.P_t = P_t;
        
        f_transition_PE('ss_autarky_q1_baseline','ss_tradeshock_q1_baseline',file_in,mat_out,t_trans,'_PE',exovars)
end




