% model_fixedcost.m.m
%
%   Compute moments for calibrated model with fixed cost here
%
%
function model_fixedcost()

%% Run model

% set path to send outputs to / xstar location
mat_out = 'fixed_cost/mat_files';
mkdir(mat_out);

% load parameters
load('common_mat_files/xstar_dM_source/xstar_fixedcost');

% run model
[~,target0,par2] = f_targets_out_fixedcost(xstart,1,1,1,1);

% compile moments
xx    = zeros(9,1);
xx(1) = target0(1);
xx(5) = target0(2);
xx(6) = target0(3);
xx(7) = target0(4);
xx(3) = target0(7);
xx(4) = 0.75*target0(8);
xx(2) = target0(5);
xx(8) = target0(6);
xx(9) = target0(9);

% export moments
moments_outcome = xx;
save([mat_out '/stats_out_fixedcost'],'moments_outcome','par2')