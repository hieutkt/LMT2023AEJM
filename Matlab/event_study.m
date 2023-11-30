function event_study(file_in)
%% event_study.m
%
%   Do simulation for event study
%%%%%%%%%%%%%

% set seed
rng('default')

% fixed R or vary R?
par_in.fixedR = 0;

% transition parameters
Tpath = 16;         % Length of panel
t_trans.T = Tpath;  % Length of panel
t_trans.fixedR = par_in.fixedR;

% number of simulated industries
N_ind = 1;


%% simulate steady model

% pre-allocate
imp_diff_data = cell(N_ind,1);
inv_data = cell(N_ind,1);
inact_data = cell(N_ind,1);
logk_data = cell(N_ind,1);
logs_data = cell(N_ind,1);
i_mpk_data = cell(N_ind,1);
id_cell = cell(N_ind,1);
year_cell = cell(N_ind,1);
imp_data_2 = cell(N_ind,1);
i_stay_data_2 = cell(N_ind,1);
k_data_2 = cell(N_ind,1);
s_data_2 = cell(N_ind,1);
i_mpk_data_2 = cell(N_ind,1);
log_mpk_i_data = cell(N_ind,1);
id_cell_2 = cell(N_ind,1);
year_cell_2 = cell(N_ind,1);
for ii=1:N_ind
    
    % simluate model
    p_stats = g_event_study(['baseline/mat_files/' file_in],t_trans);

    % save output
    inv_data{ii} = p_stats.inv_i_t;
    inact_data{ii} = p_stats.inact_i_t;
    logk_data{ii} = log(p_stats.k_i_t);
    logs_data{ii} = log(p_stats.s_i_t);
    i_mpk_data{ii} = p_stats.i_log_mpk_q_sim;
    log_mpk_i_data{ii} = p_stats.log_mpk_i_sim;
    id_cell{ii} = p_stats.id_sim;
    year_cell{ii} = p_stats.year_sim;
    
end


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
logk_data = cell2mat(reshape(logk_data,N_ind,1));
logs_data = cell2mat(reshape(logs_data,N_ind,1));
inv_rate_data = inv_data./exp(logk_data);
log_grk_rate_data = log((inv_data + (1-.105)*exp(logk_data))./exp(logk_data));
i_mpk_data = cell2mat(i_mpk_data);
log_mpk_i_data = cell2mat(log_mpk_i_data);
i_neg = (inv_rate_data<-0.01);
i_pos = (inv_rate_data>0.01);
inv_rate_data(inv_rate_data>quantile(inv_rate_data,.975)) = quantile(inv_rate_data,.975);
inv_rate_data(inv_rate_data<quantile(inv_rate_data,.025)) = quantile(inv_rate_data,.025);
    % export
p_stats_micro.micro1.imp_diff_data = imp_diff_data;
p_stats_micro.micro1.logk_data = logk_data;
p_stats_micro.micro1.logs_data = logs_data;
p_stats_micro.micro1.log_mpk_i_data = log_mpk_i_data;
p_stats_micro.micro1.inv_rate_data = inv_rate_data;
p_stats_micro.micro1.inact_data = inact_data;
p_stats_micro.micro1.ind_tag = ind_tag;
p_stats_micro.micro1.year_tag = year_tag;
p_stats_micro.micro1.id_tag = id_tag;

%% Fixed effects regression

% data for regression
xreg = [inv_rate_data log_mpk_i_data id_tag year_tag];

% sample size
Nsample = size(xreg);
Nsample = Nsample(1);

% convert to table
tbl = array2table(xreg);
tbl.Properties.VariableNames = {'inv_rate','log_mpk','id','year'};
    % list as factors
tbl.id = nominal(tbl.id);
tbl.year = nominal(tbl.year);

% declare fixed effects model: with FE
lme = fitlme(tbl,'inv_rate ~ log_mpk + (1|id)');
[beta_hat,betanames,STATS] = fixedEffects(lme);

% declare fixed effects model: without FE
lme_no_FE = fitlme(tbl,'inv_rate ~ log_mpk');
[beta_no_FE,betanames_no_FE,STATS_no_FE] = fixedEffects(lme_no_FE);

%% Plot event study

% window size
t_max = 2;

% % min years to exist in panel
% x0_trun = 2;

inv_max_store_v = zeros(3,5);
inv_min_store_v = zeros(3,5);
for x0_keep = [1:3 10]
    
    % min years to exist in panel
    x0_trun = x0_keep;
    
    unique_id = unique(p_stats_micro.micro1.id_tag);
    inv_max_store = zeros(2*t_max+1,numel(unique_id));
    inv_min_store = zeros(2*t_max+1,numel(unique_id));
    for ii=1:length(unique_id)
        ii2 = unique_id(ii);
        ix0 = find(p_stats_micro.micro1.id_tag==ii2);

        if ~isempty(ix0) && numel(ix0)>=x0_trun
            inv_rate_temp = p_stats_micro.micro1.inv_rate_data(ix0);

            %%% max
            [max_inv,ix] = max(inv_rate_temp);
            inv_max_prev = zeros(2,1);
            inv_max_next = zeros(2,1);
            for t_lag=1:t_max
                if ix>=t_lag+1
                    inv_max_prev(t_lag) = inv_rate_temp(ix-t_lag);
                else
                    inv_max_prev(t_lag) = nan;
                end
            end
            for t_next=1:t_max
                if ix<=length(inv_rate_temp)-t_next
                    inv_max_next(t_next) = inv_rate_temp(ix+1);
                else
                    inv_max_next(t_next) = nan;
                end
            end

            %%% min
            [min_inv,ix] = min(inv_rate_temp);
            inv_min_prev = zeros(2,1);
            inv_min_next = zeros(2,1);
            for t_lag=1:2
                if ix>=t_lag+1
                    inv_min_prev(t_lag) = inv_rate_temp(ix-t_lag);
                else
                    inv_min_prev(t_lag) = nan;
                end
            end
            for t_next=1:2
                if ix<=length(inv_rate_temp)-t_next
                    inv_min_next(t_next) = inv_rate_temp(ix+1);
                else
                    inv_min_next(t_next) = nan;
                end
            end

            % store the variables
            inv_max_store(:,ii) = [inv_max_prev(:) ; max_inv ; inv_max_next(:)];

            inv_min_store(:,ii) = [inv_min_prev(:) ; min_inv ; inv_min_next(:)];

        else
            % store the variables
            inv_max_store(:,ii) = nan;
            inv_min_store(:,ii) = nan;
        end
    end
    
    inv_max_store_v(x0_keep,:) = mean(inv_max_store,2,'omitnan');
    inv_min_store_v(x0_keep,:) = mean(inv_min_store,2,'omitnan');

end

%% Spike time series

% extract firm id
unique_id = unique(p_stats_micro.micro1.id_tag);

% min years to exist in panel (following Bailey-Blanco 2020)
x0_trun = 10;

% pre-allocate
t_gap = zeros(numel(unique_id),1);
t_age = zeros(numel(unique_id),1);
spike2 = zeros(numel(unique_id),1);

for ii=1:length(unique_id)
    
    % locate each firm over time
    ii2 = unique_id(ii);
    ix0 = find(p_stats_micro.micro1.id_tag==ii2);

    if ~isempty(ix0) && numel(ix0)>=x0_trun
        % some ix0 must be empty if the firm exits in 1 period

        % extract investment time series
        inv_rate_temp = p_stats_micro.micro1.inv_rate_data(ix0);

        t_spikes = find((abs(inv_rate_temp)>.49));
        if isempty(t_spikes)
            t_gap(ii) = nan;
        else
            spike2(ii) = numel(t_spikes);
            if numel(t_spikes)==1
                t_gap(ii) = length(inv_rate_temp)-t_spikes;
                t_age(ii) = t_gap(ii);
            else
                temp = zeros(numel(t_spikes)-1,1);
                for tt=1:numel(t_spikes)-1
                    temp(tt) = t_spikes(tt+1)-t_spikes(tt);
                    t_gap(ii) = mean(temp);
                end
                t_age(ii) = length(inv_rate_temp)-t_spikes(end);
            end
        end

    else
        % store the variables
        t_gap(ii) = nan;
    end
end

id_spike = spike2>=2;
t_gap_spike = t_gap(id_spike);

pop_firms_spike1 = sum(spike2>=1);
pop_firms_spike2 = sum(spike2>=2);
pop_firms = numel(t_gap);


%% Export outputs

switch file_in
    case('ss_autarky_baseline')
        save('baseline/mat_files/event_plots',...
        'inv_max_store_v','inv_min_store_v','spike2','pop_firms','pop_firms_spike1','pop_firms_spike2','t_age','t_gap','t_gap_spike',...
        'beta_hat','beta_no_FE','STATS','STATS_no_FE','Nsample');
    
        % export to stata folder
        inv_rate_975_model = inv_max_store_v(2,:)';
        event_time = [-2 -1 0 1 2]';
        T = table(event_time,inv_rate_975_model);
        writetable(T,'../stata_replication/Data/event_study_max.xlsx');
    case('ss_autarky_q1_baseline')
        save('baseline/mat_files/event_plots_q1',...
            'inv_max_store_v','inv_min_store_v','spike2','pop_firms','pop_firms_spike1','pop_firms_spike2','t_age','t_gap','t_gap_spike',...
            'beta_hat','beta_no_FE','STATS','STATS_no_FE','Nsample');
end


% print output to screen
disp('prop of firms ever having 1 spike')
disp(pop_firms_spike1/pop_firms)
disp('Average time gap till adjustment, condn on >=2 spikes')
disp(mean(t_gap_spike,'omitnan'))
disp('prop of firms having >=2 spikes')
disp(pop_firms_spike2/pop_firms)
disp('Average capital age')
disp(mean(t_age(spike2>=1),'omitnan'))
disp('Average capital age, condn on only 1 spike')
disp(mean(t_gap(spike2==1),'omitnan'))




