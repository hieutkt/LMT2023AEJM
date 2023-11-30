% model_no_eta.m
%
%   Model with convex costs, eta1,2=0
%
%   Solves the model as follows:
%   1. solve q<Q for some M that clears P
%   2. Take this M as a parameter, find P that clears q=Q
%
%   This scripts is the basis to do all the "fixed P" models with
%   transitions.
%
%   The "fixed P" version computes the counterfactual frictionless model by
%   holding fixed M solved under the baseline model, and varying the
%   market clearing price P. For the alternative (fixing P to the baseline
%   and varying M, see the "fixed M" version.
%
%   It solves the model as follows:
%       1. solve q<Q for some M that clears P (calls f_PE_ss_autarky)
%           a. Given M, solve trade shock model by varying P (calls f_PE_ss_trade_shock)
%           b. Taking the two end points as given, solve for the full transition path (calls f_transition)
%       2. solve q=Q by taking in M as given, and clearing the market using P (calls f_ss_autarky)
%           a. Given M, solve trade shock model by varying P (calls f_PE_ss_trade_shock)
%           b. Taking the two end points as given, solve for the full transition path (calls f_transition)
%
%   Note: Use switch to switch between github and dropbox
%
function model_no_eta(model)

%% Housekeeping
clc
close all

% set path to send outputs to / xstar location
mat_out = 'no_eta/mat_files';
mkdir(mat_out);

% to save output or not?
par_in.save = 1;

% calibrating or actually solving model?
par_in.cali = 0;

% fixed R (PE) or vary R (GE)?
par_in.fixedR = 0;

% baseline fixed parameters 
par_in.rho_fstay = 0;
par_in.Me = 5;
switch model
    case 'baseline'
        par_in.Ns = 31;
    case 'frictionless'
        par_in.Ns = 31*3;
end
par_in.Nk = 1000;
% par_in.Nk = 10000;
par_in.new_kmax = 8.0526;
par_in.new_Pnow = 1/3; % fixed P (Clementi+Palazzo calibration)
par_in.Q = par_in.new_Pnow;
par_in.npow = 1;
par_in.kapp = 1;
par_in.s_sh = 1;
par_in.mu_fstay_a = 0;
par_in.mu_fenter_a = 0;

% load calibrated baseline parameters
load('common_mat_files/xstar_dM_source/xstar_no_eta.mat', 'xout')
x0 = xout;
xstar(1)            = x0(1);
xstar(2)            = x0(2);
par_in.c1           = x0(3);
par_in.eta(1)       = 0;
par_in.eta(2)       = 0;
par_in.mu_fstay_b   = (par_in.mu_fstay_a+x0(4));
par_in.mu_fenter_b  = (par_in.mu_fenter_a+x0(4));
par_in.rho_s        = x0(5);
par_in.sigma_us     = x0(6);
    % impose identical convex cost for intensive / extensive margins
par_in.c2 = par_in.c1;
    % fixed cost (==0 in baseline)
par_in.C_f = 0;
    % lam and zet
switch model
    
    case 'baseline'
        % q<Q , zet>0
        par_in.lam = xstar(1);
        par_in.zet = xstar(2);
        
        disp('Solving model with frictions')
        
    case 'frictionless'
        % q=Q , zet=0
        par_in.lam = 0;
        par_in.zet = 0;
        
        disp('Solving model without frictions')
        
end

%% Run the steady states and transitions

%%% solve closed economy steady-state
switch model
    case 'baseline'
        
        closed_moments = f_PE_ss_autarky_convex(par_in,mat_out,'_no_eta');
        
        % save measure of entrants
        load([mat_out '/ss_autarky_no_eta'], 'Me_new');
        Mnow = Me_new;
        save([mat_out '/Mnow'], 'Mnow')
        
    case{'frictionless'}
        % load the correct measure
        try
            load([mat_out '/Mnow'],'Mnow')
        catch
            disp('Need to solve q<Q first!')
        end
        par_in.Me = Mnow;
        closed_moments = f_ss_autarky(par_in,mat_out,'_no_eta');

        clear Mnow
end

%%% solve open economy steady-state
% add parameters for open economy and solve
switch model
    case 'baseline'
        load([mat_out '/Pss_autarky_no_eta'], 'Pnow')

    case 'frictionless'
        load([mat_out '/Pss_autarky_q1_no_eta'], 'Pnow')     
end
    % load calibrated mass of entrants from autarky
load([mat_out '/Mnow'],'Mnow')
    % load calibrated mass of foreign firms
load('common_mat_files/xstar_dM_source/dMf_star.mat', 'dMf_star')
par_in.dM           = dMf_star;     % mass of foreign firms
par_in.Pf           = 0.5*Pnow;     % price of foreign goods
par_in.P_autarky    = Pnow;
par_in.Me           = Mnow;         % update with measure from closed economy
    % impose kmax at value set in autarky solution
switch model
    case 'baseline'
        
        load([mat_out '/kmax_autarky_no_eta'], 'kmax_autarky')
        
    case 'frictionless'
        
        load([mat_out '/kmax_autarky_q1_no_eta'], 'kmax_autarky')
end
par_in.kmax_autarky = kmax_autarky;
clear Pnow kmax_autarky Mnow
    % solve model
switch model
    case{'baseline'}
        opened_moments = f_PE_ss_trade_shock_convex(par_in,mat_out,'_no_eta');
        
    case{'frictionless'}
        opened_moments = f_PE_ss_trade_shock(par_in,mat_out,'_no_eta');
end

%%% solve transition path
t_trans.Pf = 0;
t_trans.T = 1;
t_trans.Tmax = 40;
t_trans.fixedR = par_in.fixedR;
t_trans.Q = par_in.Q;
switch model
    case 'baseline'
        f_transition_mkt_convex('ss_autarky_no_eta','ss_tradeshock_no_eta',mat_out,t_trans,'_no_eta')
        
    case 'frictionless'
        f_transition_mkt('ss_autarky_q1_no_eta','ss_tradeshock_q1_no_eta',mat_out,t_trans,'_no_eta')
end




