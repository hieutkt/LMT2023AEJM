% model_recreate_transition.m
%
%   Feed in path of prices from GE (or PE) to recreate transition
%   statistics and compute decomposition exercises
%
function model_recreate_transition(model,calibration_version)

%% Housekeeping
clc
close all

% set path to send outputs to / xstar location
mat_out = 'Recreated_workspace/mat_files';
mkdir(mat_out);

% fixed R (PE) or vary R (GE)?
par_in.fixedR = 0;

par_in.new_Pnow = 1/3; % fixed P (Clementi+Palazzo calibration)
par_in.Q = par_in.new_Pnow;

%% Run the steady states and transitions

%%% solve transition path
t_trans.Pf = 0;
t_trans.T = 1;
t_trans.Tmax = 40;
t_trans.fixedR = par_in.fixedR;
t_trans.Q = par_in.Q;

switch calibration_version
    case('baseline calibration')
        
        % path to solution
        file_in = 'mat_files_compiled/baseline/mat_files/';
        
        % switch between baseline and frictionless
        switch model
            case 'baseline'
                load('mat_files_compiled/baseline/mat_files/transition_baseline.mat','P_t')
                exovars.P_t = P_t;

                f_recreate_transition_PE_convex('ss_autarky_baseline','ss_tradeshock_baseline',file_in,mat_out,t_trans,'_baseline',exovars)

            case 'frictionless'
                load('mat_files_compiled/baseline/mat_files/transition_q1_baseline.mat','P_t')
                exovars.P_t = P_t;

                f_recreate_transition_PE('ss_autarky_q1_baseline','ss_tradeshock_q1_baseline',file_in,mat_out,t_trans,'_baseline_q1',exovars)
        end
        
    case('no eta calibration')
        
        % path to solution
        file_in = 'mat_files_compiled/no_eta/mat_files/';
        
        % switch between baseline and frictionless
        switch model
            case 'baseline'
                load('mat_files_compiled/no_eta/mat_files/transition_no_eta.mat','P_t')
                exovars.P_t = P_t;

                f_recreate_transition_PE_convex('ss_autarky_no_eta','ss_tradeshock_no_eta',file_in,mat_out,t_trans,'_no_eta',exovars)

            case 'frictionless'
                load('mat_files_compiled/no_eta/mat_files/transition_q1_no_eta.mat','P_t')
                exovars.P_t = P_t;

                f_recreate_transition_PE('ss_autarky_q1_no_eta','ss_tradeshock_q1_no_eta',file_in,mat_out,t_trans,'_no_eta_q1',exovars)
        end
        
end



