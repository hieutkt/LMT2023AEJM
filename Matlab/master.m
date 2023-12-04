
% master.m
%
%   Master file for running all models
%
%   For each model, run the "baseline" before "frictionless"
%

disp("[master.m] Starts replicating, turning on timer...")
%%%%%%%%%%%%%

clear
clc
close all
delete(gcp('nocreate'));

% set path to functions
addpath('functions')
addpath('utility functions');

% removed mat_files_compiled folder to avoid file conflict
try
    rmdir('mat_files_compiled','s')
    disp('    mat_files_compiled exist, removing folder to avoid path conflict')
catch
    disp('    mat_files_compiled does not already exist')
end
mkdir('mat_files_compiled')

disp('.')
pause(0.5)
disp('.')
pause(0.5)
disp('.')
pause(0.5)

% baseline model
disp("[master.m] Solving the autarky model; Entering model_baseline.m ...")
model_baseline('baseline')

disp("[master.m] Solving the frictionless model; Entering model_baseline.m ...")
model_baseline('frictionless')

disp("[master.m] Solving the covid model; Entering model_baseline.m ...")
model_baseline('covid')


disp("[master.m] Solving the covid model; Entering model_baseline.m ...")
model_baseline('covid_frictionless')


disp("[master.m] Extract the outputs...")
extract_outputs('baseline/mat_files', ...
                'ss_autarky_baseline','ss_tradeshock_baseline', 'transition_baseline', ...
                'ss_autarky_q1_baseline','ss_tradeshock_q1_baseline', 'transition_q1_baseline', ...
                'uniform')

% extract_outputs('baseline/mat_files', ...
%                 'ss_tradeshock_baseline', 'ss_covid_baseline', 'transition_baseline', ...
%                 'ss_tradeshock_q1_baseline', 'ss_covid_q1_baseline', 'transition_q1_baseline', ...
%                 'uniform')


%     % run event study + i/k vs mrpk elasticity
% disp("[master.m] run event study + i/k vs mrpk elasticity")
% event_study('ss_autarky_baseline')
% event_study('ss_autarky_q1_baseline')

    % move outputs to combined folder
disp("[master.m] move outputs to combined folder")
i_success = movefile('baseline','mat_files_compiled');

%     % decomposition exercise: baseline calibration
% disp("[master.m] decomposition exercise: baseline calibration")
% model_recreate_transition('baseline','baseline calibration')
% model_recreate_transition('frictionless','baseline calibration')
%     % multi-industry trade shock (serial)
% disp("[master.m] multi-industry trade shock (serial)")
% model_trade_shock_baseline_serial('mat_files_compiled/baseline/mat_files')
% model_trade_shock_frictionless_serial('mat_files_compiled/baseline/mat_files')
%     % compile event study mat files with multi-industry trade shock
% disp("[master.m] compile event study mat files with multi-industry trade shock")
% i_success = movefile('mat_files_compiled/baseline/mat_files/event_plots.mat','mat_files_compiled/baseline/mat_files/results');
% i_success = movefile('mat_files_compiled/baseline/mat_files/event_plots_q1.mat','mat_files_compiled/baseline/mat_files/results');

% % PE model
% model_baseline_PE('baseline')
% model_baseline_PE('frictionless')
% extract_outputs('baseline_PE/mat_files','../../mat_files_compiled/baseline/mat_files/ss_autarky_baseline','../../mat_files_compiled/baseline/mat_files/ss_tradeshock_baseline',...
%     'transition_PE','../../mat_files_compiled/baseline/mat_files/ss_autarky_q1_baseline',...
%     '../../mat_files_compiled/baseline/mat_files/ss_tradeshock_q1_baseline','transition_q1_PE','uniform')
%     % move outputs to combined folder
% i_success = movefile('baseline_PE','mat_files_compiled');

% % eta_1=eta_2=0 model
% model_no_eta('baseline')
% model_no_eta('frictionless')
% extract_outputs('no_eta/mat_files','ss_autarky_no_eta','ss_tradeshock_no_eta','transition_no_eta','ss_autarky_q1_no_eta','ss_tradeshock_q1_no_eta','transition_q1_no_eta','uniform')
%     % move outputs to combined folder
% i_success = movefile('no_eta','mat_files_compiled');
%     % decomposition exercise: no eta
% model_recreate_transition('baseline','no eta calibration')
% model_recreate_transition('frictionless','no eta calibration')

% % high sigma model
% model_sigstar('baseline')
% model_sigstar('frictionless')
% extract_outputs('high_sigma/mat_files','ss_autarky_high_sigma','ss_tradeshock_high_sigma','transition_high_sigma','ss_autarky_q1_high_sigma','ss_tradeshock_q1_high_sigma','transition_q1_high_sigma','uniform')
%     % move outputs to combined folder
% i_success = movefile('high_sigma','mat_files_compiled');

% % no convex cost model
% model_no_convex('baseline')
% model_no_convex('frictionless')
% extract_outputs('no_convex/mat_files','ss_autarky_no_convex','ss_tradeshock_no_convex','transition_no_convex','ss_autarky_q1_no_convex','ss_tradeshock_q1_no_convex','transition_q1_no_convex','uniform')
%     % move outputs to combined folder
% i_success = movefile('no_convex','mat_files_compiled');

% % fixed cost model
%     % no inputs needed
% model_fixedcost
%     % move outputs to combined folder
% i_success = movefile('fixed_cost','mat_files_compiled');

% move recreated workspace to combined folder
% disp("[master.m] move recreated workspace to combined folder")
% i_success = movefile('Recreated_workspace','mat_files_compiled');

% Plot / construct all figures and tables
disp("[master.m] Plot / construct all figures and tables")
CHY_gen_figures_tables

% generate macros used in table
% disp("[master.m] generate macros used in table")
% gen_macros

% export selected simulation to stata
% disp("[master.m] export selected simulation to stata")
% tostata
