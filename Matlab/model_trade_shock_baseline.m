%% pame_stats_script.m
%
%   Generate regression stats with import penetration shocks for CONVEX
%   model
%
%   Code automatically loads initial / final / transition workspaces
%
%   v1: "GE" version.
%   v2: "Andrea" version
%%%%%%%%%%%%%

function model_trade_shock_baseline(mat_out)
clc
close all


% fix seed
rng('default')

% directory to export results
mkdir([mat_out '/results'])

% fixed R or vary R?
par_in.fixedR = 1;

% transition parameters
t_trans.Pf = 0;
t_trans.T = 1;
t_trans.Tmax = 40;
t_trans.fixedR = par_in.fixedR;
    % years to extract
Tpath = 16;
    % number of simulated industries
N_ind = 8*4;


% "persistent shocks"
imp_shocks_draw = zeros(16,N_ind);
t_end = 2;
delta_p = 0.;
for tt = 1:Tpath
    
    if tt==1
        % initial impulse
%         imp_shocks_draw(1,:) = 0.1*randn(1,N_ind);
        imp_shocks_draw(1,:) = 0.065*(rand(1,N_ind)-0.5)/0.5;
        
    elseif tt>1 && tt<=t_end
        
        % dependence until t_end
        imp_shocks_draw(tt,:) = delta_p^(tt-1)*imp_shocks_draw(1,:);
    else
        
        % fade out
        imp_shocks_draw(tt,:) = 0;
    end
    
end

%% v1: GE 

clear bhat_out stats_agg stats_micro

t_trans.andrea = 0;

import_pen_t_realized = zeros(Tpath-1,N_ind);
    % agg stats
sd_mpk_t = zeros(Tpath-1,N_ind);
TFPQ_t = zeros(Tpath-1,N_ind);
inact_t = zeros(Tpath-1,N_ind);
imp_t = zeros(Tpath-1,N_ind);
inv_data = cell(N_ind,1);
inact_data = cell(N_ind,1);
ipos_t = zeros(Tpath-1,N_ind);
ineg_t = zeros(Tpath-1,N_ind);
avg_ik_t = zeros(Tpath-1,N_ind);
agg_ik_t = zeros(Tpath-1,N_ind);
P_sim_t = zeros(Tpath-1,N_ind);
exit_rate_t = zeros(Tpath-1,N_ind);
    % firm level
imp_diff_data = cell(N_ind,1);
imp_pen_data = cell(N_ind,1);
logk_data = cell(N_ind,1);
logs_data = cell(N_ind,1);
i_mpk_data = cell(N_ind,1);
log_mpk_data = cell(N_ind,1);
id_cell = cell(N_ind,1);
year_cell = cell(N_ind,1);
imp_data_2 = cell(N_ind,1);
i_stay_data_2 = cell(N_ind,1);
k_data_2 = cell(N_ind,1);
s_data_2 = cell(N_ind,1);
i_mpk_data_2 = cell(N_ind,1);
id_cell_2 = cell(N_ind,1);
year_cell_2 = cell(N_ind,1);
logk_t1 = cell(N_ind,1);
logs_t1 = cell(N_ind,1);
inact_t1 = cell(N_ind,1);
inv_rate_t1 = cell(N_ind,1);
i_stay_t1 = cell(N_ind,1);
i_pos_t1 = cell(N_ind,1);
i_neg_t1 = cell(N_ind,1);
logk_all_t1 = cell(N_ind,1);
logs_all_t1 = cell(N_ind,1);
istay_all_t1 = cell(N_ind,1);
parfor ii=1:N_ind
% parfor ii=1:N_ind
    
    % generate stats
    out_stats = h_transition_shocks_convex('ss_autarky_baseline','ss_tradeshock_baseline','transition_baseline',imp_shocks_draw(:,ii),mat_out,t_trans);
    
    % aggregate industry data
    import_pen_t_realized(:,ii) = out_stats.import_pen_t;
    sd_mpk_t(:,ii) = out_stats.sd_mpk_t;
    TFPQ_t(:,ii) = out_stats.TFPQ_t;
    inact_t(:,ii) = out_stats.inact_t;
    ipos_t(:,ii) = out_stats.ipos;
    ineg_t(:,ii) = out_stats.ineg;
    avg_ik_t(:,ii) = out_stats.avg_ik;
    agg_ik_t(:,ii) = out_stats.agg_ik;
    P_sim_t(:,ii) = out_stats.P_t;
    exit_rate_t(:,ii) = out_stats.exit_rate_t;
    
    % micro full time series, continuers
    inv_data{ii} = out_stats.inv_i_t;
    inact_data{ii} = out_stats.inact_i_t;
    imp_diff_data{ii} = out_stats.imp_pen_diff_sample;
    imp_pen_data{ii} = out_stats.imp_pen_sample;
    logk_data{ii} = log(out_stats.k_i_t);
    logs_data{ii} = log(out_stats.s_i_t);
    i_mpk_data{ii} = out_stats.i_log_mpk_q_sim;
    log_mpk_data{ii} = out_stats.log_mpk_i_sim;
        % id tags
    id_cell{ii} = out_stats.id_sim;
    year_cell{ii} = out_stats.year_sim;
    
    % micro full time series, all firms
    imp_data_2{ii} = out_stats.imp_pen_diff_sample_2;
    imp_pen_data_2{ii} = out_stats.imp_pen_sample_2;
    i_stay_data_2{ii} = out_stats.i_stay_vec;
    k_data_2{ii} = out_stats.k_active_sim;
    s_data_2{ii} = out_stats.s_active_sim;  
    i_mpk_data_2{ii} = out_stats.i_log_mpk_all_q_sim;
        % id tags
    year_cell_2{ii} = out_stats.year_allfirms_sim;
    id_cell_2{ii} = out_stats.id_allfirms_sim;

    % t = 1 data, continuers only
    logk_t1{ii} = log(out_stats.kstate_vec_t1);
    logs_t1{ii} = log(out_stats.s_vec_t1);
    inact_t1{ii} = out_stats.inact_i_t1;
    inv_rate_t1{ii} = out_stats.inv_i_rate_t1;
    i_stay_t1{ii} = out_stats.i_stay_sim_t1;
    i_pos_t1{ii} = out_stats.i_pos_t1;
    i_neg_t1{ii} = out_stats.i_neg_t1;
    
    % t = 1 data, all firms
    logk_all_t1{ii} = out_stats.k_active_sim_t1;
    logs_all_t1{ii} = out_stats.s_active_sim_t1;
    istay_all_t1{ii} = out_stats.i_stay_sim_t1;
    
end

% construct agg equivalent of imp diff
load([mat_out '/transition_baseline'],'P_Cd_t','P_Cf_t')
import_pen_trend = 1-(P_Cd_t(1:Tpath-1) ./ (P_Cd_t(1:Tpath-1) + P_Cf_t(1:Tpath-1)));
import_pen_diff_t = import_pen_t_realized - import_pen_trend;
import_pen_diff_t1 = import_pen_diff_t(1,:);
import_pen_t1_realized = import_pen_t_realized(1,:);

% agg stats
stats_agg.TFPQ_t = TFPQ_t;
stats_agg.sd_mpk_t = sd_mpk_t;
stats_agg.inact_t = inact_t;
stats_agg.ipos_t = ipos_t;
stats_agg.ineg_t = ineg_t;
stats_agg.import_pen_t_realized = import_pen_t_realized;
stats_agg.import_pen_diff_t = import_pen_diff_t;
stats_agg.imp_shocks_draw = imp_shocks_draw;
stats_agg.avg_ik_t = avg_ik_t;
stats_agg.agg_ik_t = agg_ik_t;
stats_agg.P_sim_t = P_sim_t;
stats_agg.exit_rate_t = exit_rate_t;

% full time series: continuers only
    % industry tags
temp_cell = cell(N_ind,1);
for ii=1:N_ind
   temp_cell{ii} = repmat(ii,size(inv_data{ii})); 
end
ind_tag = cell2mat(temp_cell);
clear temp_cell
    % year tags
year_tag = cell2mat(reshape(year_cell,N_ind,1));
    % id tags
id_tag = cell2mat(reshape(id_cell,N_ind,1));
    % convert data to vector
inv_data = cell2mat(inv_data);
inact_data = cell2mat(inact_data);
imp_diff_data = cell2mat(imp_diff_data);
imp_pen_data = cell2mat(imp_pen_data);
logk_data = cell2mat(reshape(logk_data,N_ind,1));
logs_data = cell2mat(reshape(logs_data,N_ind,1));
inv_rate_data = inv_data./exp(logk_data);
log_grk_rate_data = log((inv_data + (1-.105)*exp(logk_data))./exp(logk_data));
i_mpk_data = cell2mat(i_mpk_data);
log_mpk_data = cell2mat(log_mpk_data);
i_neg = (inv_rate_data<-0.01);
i_pos = (inv_rate_data>0.01);
    % export
stats_micro.micro1.imp_diff_data = imp_diff_data;    
stats_micro.micro1.imp_pen_data = imp_pen_data;
stats_micro.micro1.logk_data = logk_data;
stats_micro.micro1.logs_data = logs_data;
stats_micro.micro1.inv_rate_data = inv_rate_data;
stats_micro.micro1.inact_data = inact_data;
stats_micro.micro1.log_mpk_data = log_mpk_data;
stats_micro.micro1.ind_tag = ind_tag;
stats_micro.micro1.year_tag = year_tag;
stats_micro.micro1.id_tag = id_tag;

% full time series: all firms
    % industry tags
temp_cell = cell(N_ind,1);
for ii=1:N_ind
   temp_cell{ii} = repmat(ii,size(imp_data_2{ii})); 
end
ind_tag_2 = cell2mat(temp_cell);
clear temp_cell
    % year tags
year_tag_2 = cell2mat(reshape(year_cell_2,N_ind,1));
    % id tags
id_tag_2 = cell2mat(reshape(id_cell_2,N_ind,1));
    % convert data to vector
imp_data_2 = cell2mat(reshape(imp_data_2,N_ind,1));
imp_pen_data_2 = cell2mat(reshape(imp_pen_data_2,N_ind,1));
k_data_2 = cell2mat(reshape(k_data_2,N_ind,1));
s_data_2 = cell2mat(reshape(s_data_2,N_ind,1));
i_stay_data_2  = cell2mat(reshape(i_stay_data_2,N_ind,1));
i_mpk_data_2 = cell2mat(i_mpk_data_2);
    % export
stats_micro.micro2.imp_diff_data = imp_data_2;
stats_micro.micro2.imp_pen_data_2 = imp_pen_data_2;
stats_micro.micro2.logk_data = log(k_data_2);
stats_micro.micro2.logs_data = log(s_data_2);
stats_micro.micro2.i_stay_data = i_stay_data_2;
stats_micro.micro2.i_mpk_data = i_mpk_data_2;
stats_micro.micro2.ind_tag = ind_tag_2;
stats_micro.micro2.year_tag = year_tag_2;
stats_micro.micro2.id_tag = id_tag_2;
% % % inv_vec = zeros(size(imp_data_2));
% % % inv_vec(i_stay_data_2) = log_grk_rate_data;
% % % inv_vec(~i_stay_data_2) = -1e10;
% % % stats_micro.micro2.inv_rate_data = inv_vec;

% for t=1
    % continuers only
stats_micro.t_1_cont.logk_t1 = cell2mat(logk_t1);
stats_micro.t_1_cont.logs_t1 = cell2mat(logs_t1);
stats_micro.t_1_cont.inact_t1 = cell2mat(inact_t1); % inaction == 1
stats_micro.t_1_cont.inv_rate_t1 = cell2mat(inv_rate_t1); 
stats_micro.t_1_cont.i_stay_t1 = cell2mat(i_stay_t1); % stay == 1
stats_micro.t_1_cont.i_pos_t1 = cell2mat(i_pos_t1); % pos == 1
stats_micro.t_1_cont.i_neg_t1 = cell2mat(i_neg_t1); % neg == 1
temp_cell = cell(N_ind,1);
for ii=1:N_ind
   temp_cell{ii} = repmat(import_pen_diff_t1(ii),size(logk_t1{ii})); 
end
stats_micro.t_1_cont.import_pen_diff_t1 = cell2mat(temp_cell);
clear temp_cell
    % all firms
stats_micro.t_1_all.logk_all_t1 = log(cell2mat(logk_all_t1));
stats_micro.t_1_all.logs_all_t1 = log(cell2mat(logs_all_t1));
stats_micro.t_1_all.istay_all_t1 = cell2mat(istay_all_t1);
% stats_micro.t_1_all.import_pen_diff_t1 = import_pen_diff_t1;
temp_cell = cell(N_ind,1);
temp_cell2 = cell(N_ind,1);
for ii=1:N_ind
   temp_cell{ii} = repmat(import_pen_diff_t1(ii),size(logk_all_t1{ii})); 
   temp_cell2{ii} = repmat(import_pen_t1_realized(ii),size(logk_all_t1{ii})); 
end
stats_micro.t_1_all.import_pen_diff_t1 = cell2mat(temp_cell);
stats_micro.t_1_all.import_pen_t1 = cell2mat(temp_cell2);
clear temp_cell

save([mat_out '/results/stats_data'],'stats_micro','stats_agg','imp_shocks_draw');

%%%%%%%% Results 

% aggregates
bhat_out.tfp = regress(log(TFPQ_t(:)),[ones(numel(import_pen_diff_t),1) import_pen_diff_t(:)]);
bhat_out.sdpmk = regress(sd_mpk_t(:),[ones(numel(import_pen_diff_t),1) import_pen_diff_t(:)]);
bhat_out.inact = regress(inact_t(:),[ones(numel(import_pen_diff_t),1) import_pen_diff_t(:)]);

% results with only continuers
bhat_out.invrate_micro = regress(inv_rate_data(:),[ones(numel(imp_diff_data),1) imp_diff_data(:) logk_data(:) logs_data(:)]);
bhat_out.inact_micro = regress(inact_data(:),[ones(numel(imp_diff_data),1) imp_diff_data(:) logk_data(:) logs_data(:)]);
bhat_out.i_neg_micro = regress(i_neg(:),[ones(numel(imp_diff_data),1) imp_diff_data(:) logk_data(:) logs_data(:)]);
bhat_out.i_pos_micro = regress(i_pos(:),[ones(numel(imp_diff_data),1) imp_diff_data(:) logk_data(:) logs_data(:)]);
[bhat_out.invrate_x_mpk_micro,ci_out.invrate_x_mpk_micro] = regress(log_grk_rate_data(:),[ones(numel(imp_diff_data),1) ...
    (i_mpk_data(:)==2) (i_mpk_data(:)==3) ...
    imp_diff_data(:).*(i_mpk_data(:)==1) ...
    imp_diff_data(:).*(i_mpk_data(:)==2) imp_diff_data(:).*(i_mpk_data(:)==3) ...
    logk_data(:) logs_data(:)]);

% implied persistence of import penetration shock
import_pen_diff_t_now = stats_agg.import_pen_diff_t(2:end,:);
import_pen_diff_t_lag = stats_agg.import_pen_diff_t(1:end-1,:);
bhat_out.AR1 = regress(import_pen_diff_t_now(:),[ones(size(import_pen_diff_t_lag(:))) import_pen_diff_t_lag(:)]);

% exit regression, full sample
bhat_out.cont_reg = glmfit([imp_data_2(:) log(k_data_2(:)) log(s_data_2(:)) ...
    imp_data_2(:).*log(s_data_2(:)) ...
     imp_data_2(:).*log(k_data_2(:)) ...
    ],i_stay_data_2(:),'binomial','link','probit');

% exit regression, full sample
bhat_out.cont_reg = glmfit([imp_data_2(:) log(k_data_2(:)) log(s_data_2(:)) ...
    imp_data_2(:).*log(s_data_2(:)) ...
     imp_data_2(:).*log(k_data_2(:)) ...
    ],i_stay_data_2(:),'binomial','link','probit');

% one-period exit regression
bhat_out.bhat_exit_t1 = glmfit([stats_micro.t_1_all.import_pen_t1(:) (stats_micro.t_1_all.logk_all_t1(:)) (stats_micro.t_1_all.logs_all_t1(:)) ...
    ],stats_micro.t_1_all.istay_all_t1(:),'binomial','link','probit');

bhat_out.bhat_exit_logs_t1 = glmfit([log(stats_micro.t_1_all.import_pen_t1(:)) (stats_micro.t_1_all.logk_all_t1(:)) (stats_micro.t_1_all.logs_all_t1(:)) ...
    ],stats_micro.t_1_all.istay_all_t1(:),'binomial','link','probit');

% export regression results
save([mat_out '/results/regs'],'bhat_out','ci_out')

