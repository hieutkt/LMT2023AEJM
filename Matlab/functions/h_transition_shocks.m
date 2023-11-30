function out_stats = h_transition_shocks(init_ss,final_ss,transition_file,p_shocks,file_out,trans_x,file_suffix)
%

% path to workspaces
if nargin==2
    file_in = [];
elseif nargin>=3
    file_in = [file_out '/'];
end

% load transition file first
load([file_in transition_file], 'P_t', 'Pf_t','dM_t')
T = min(16,length(P_t));
P_t = P_t(1:T);
Pf_t = Pf_t(1:T);
dM_t = dM_t(1:T);

% storage
vf_t = cell(T,1);
pol_t = cell(T,1);
gr_t = cell(T,1);

% load initial steady state
load([file_in init_ss], 'Px', 'par2', 'gr', 'N', 'Pnow', 'old_C','pol')
% load([file_in init_ss], 'Pnow', 'old_C'); % debug
% load([file_in final_ss], 'Px', 'par2', 'gr', 'N'); % debug
g0_autarky = squeeze(Px.g_ss);
P_autarky = Pnow;
C_autarky = old_C;
par_t0 = par2;
kmax_autarky = gr.k_grid(end);
par_t0.kstar = kmax_autarky;
par_t0.P0 = P_autarky;
par_t0.C0 = C_autarky;
gr0.k_grid = gr.k_grid;
pol_0 = pol;
% when debugging
% try
%     load([file_in init_ss], 'Pss_notrade')
%     Q_autarky = Pss_notrade;
%     clear Pss_notrade
%     
% catch
    Q_autarky = par2.Q;
    
% end
clear Px Pnow par gr new_C pol

% load final steady state
load([file_in final_ss], 'initval', 'pol', 'par2', 'Pf', 'dM', 'Pnow', 'old_C')
% load([file_in init_ss], 'initval', 'pol', 'par2', 'Pf', 'dM', 'Pnow', 'old_C'); % debug
vf_t{end} = initval.Ev;
pol_t{end} = pol;
P_trade = Pnow;
C_trade = old_C;
par_tT = par2;
Pf_trade = Pf;
dM_trade = dM;
clear initval pol Pnow par new_C Pf dM


% trans_x contains switch for fixedR since par_in is not passed in here
par_t0.fixedR = trans_x.fixedR;

% save trend prices
P_t_input = P_t;

% allocate shocked prices
P_t = zeros(size(p_shocks));
    % initial impulse
% P_t(1,:) = P_t_input(1).*exp(p_shocks(1,:));
P_t(1,:) = P_t_input(1).*(1+p_shocks(1,:));
for t=2:length(P_t_input)
    xt = P_t(1,:)./P_t_input(1);
    if t==2
        % impose immediate mean reversion
        P_t(t,:) = P_t_input(2)*xt;
    else
        P_t(t,:) = P_t_input(t);
    end
end

% Back out "GE" interest rates
    % Andrea's version where SDF in labor units move around (wrong, don't use)
% R_t = (P_t_input(2:end)./P_t_input(1:end-1)) .* (P_t(1:end-1)./P_t(2:end));
R_t = P_t(1:end-1)./P_t(2:end);
% R_t = ones(size(R_t));
if trans_x.andrea==0
    % my version where SDF in labor units is fixed
    R_t = ones(size(R_t));
end


%% Solve for policy functions 

for t=T-1:-1:1

    % always reset to t0 par
    par2 = par_t0;

    % next period value function
    initval.Ev = vf_t{t+1};
    % next period price
    par2.P_pr = P_t(t+1);
    % this period price
    par2.P = P_t(t);
    % last period price
    par2.P_last = P_t(max(1,t-1));
    % fix price of capital
    par2.Q = Q_autarky;
    % old GE price
    par2.P_old = P_t_input(t);
    par2.P_pr_old = P_t_input(t+1);
    % fix consumption aggregate
    par2.C0 = 1/(par2.P_old*par2.chi);
    par2.C0_pr = 1/(par2.P_pr_old*par2.chi);
    % foreign price
    Pf = Pf_t(t);
    % foreign mass of firms
    dM = dM_t(t);
    % time / t0 k_grid
    if par2.lam==0 && par2.zet==0
        par2.t = t;
        par2.k_grid = gr0.k_grid;
    end
    %%%%%%%%%%
    %  need to edit for set_grid vs policy part (for q=Q, the interest rate
    %  denoted in C units in set_grid is needed)
    % Interest rate
    par2.R = R_t(t);
    if t==1
        par2.R_last = 1;
    else
        par2.R_last = R_t(t-1);
    end

    %%% set the grid
    gr = mainfuncMP_fixed_C(par2,N,[],[],2);

    %%% solve the individual's problem by backward induction
    output  = mainfuncMP_fixed_C(par2,N,gr,[],5,initval);
    pol     = output.pol;
    v       = output.vf;
    clear output

    gr.i_inact_xgrid = [];
    gr.k_inact_xgrid = [];
    gr.i_inact_xgrid = [];
    gr.s_xgrid0 = [];
    gr.tau_xgrid = [];
    gr.x_xgrid = [];
    gr.s_xgrid = [];
    gr.py_i = [];
    gr.y_tau = [];
    if t>1
        gr.k_xgrid = [];
    end

    vf_t{t} = v;
    pol_t{t} = pol;
    gr_t{t} = gr;

end


%% Iterate forward, save stats

C_t = zeros(T-1,1);
Cd_t = zeros(T-1,1);
Cf_t = zeros(T-1,1);
P_Cd_t = zeros(T-1,1);
P_Cf_t = zeros(T-1,1);
Kd_t = zeros(T-1,1);
Ld_t = zeros(T-1,1);
L_t = zeros(T-1,1);
TFPR_t = zeros(T-1,1);
TFPQ_t = zeros(T-1,1);
Mentrants_t = zeros(T-1,1);
Mactive_t = zeros(T-1,1);
Mexit_t = zeros(T-1,1);
M_foreign = zeros(T-1,1);

% Aggregates
IK_t = zeros(T-1,1);
I_Gross_t = zeros(T-1,1);
I_Gross_plus_t = zeros(T-1,1);
I_Gross_neg_t = zeros(T-1,1);
I_net_t = zeros(T-1,1);

% Expected values
I_Gross_intensive_t = zeros(T-1,1);
I_Gross_plus_intensive_t = zeros(T-1,1);
I_Gross_neg_intensive_t = zeros(T-1,1);
I_Gross_frac_0_intensive_t = zeros(T-1,1);
I_Gross_frac_pos_intensive_t = zeros(T-1,1);
I_Gross_frac_neg_intensive_t = zeros(T-1,1);
I_Gross_frac_intensive_neg_t = zeros(T-1,1);

operating_cost_t = zeros(T-1,1);
entry_cost_t = zeros(T-1,1);

avg_raw_s_t = zeros(T-1,1);
avg_wgt_s_t = zeros(T-1,1);
avg_wgt_y_s_t = zeros(T-1,1);
avg_wgt_k_s_t = zeros(T-1,1);

avg_raw_s_entrants_t = zeros(T-1,1);
avg_raw_s_exiters_t = zeros(T-1,1);
avg_raw_s_cont_t = zeros(T-1,1);

avg_s_med = zeros(T-1,1);
avg_s_med_2 = zeros(T-1,1);

mean_mpk_t = zeros(T-1,1);
sd_mpk_t = zeros(T-1,1);
sd_log_s_t = zeros(T-1,1);

avg_k_exit_t = zeros(T-1,1);
avg_k_cont_t = zeros(T-1,1);
avg_k_enter_t = zeros(T-1,1);
avg_frac_neg_kpr_t = zeros(T-1,1);
    
g0 = g0_autarky;
clear initval
for t=1:T-1

    % always reset to t0 par
    par2 = par_t0;
    % next period price
    par2.P_pr = P_t(t+1);
    % this period price
    par2.P = P_t(t);
    % last period price
    par2.P_last = P_t(max(1,t-1));
    % fix price of capital
    par2.Q = Q_autarky;
    % old GE price
    par2.P_old = P_t_input(t);
    par2.P_pr_old = P_t_input(t+1);
    % fix consumption aggregate
    par2.C0 = 1/(par2.P_old*par2.chi);
    par2.C0_pr = 1/(par2.P_pr_old*par2.chi);
    % foreign price
    Pf = Pf_t(t);
    % foreign mass of firms
    dM = dM_t(t);
    % time / t0 k_grid
    if par2.lam==0 && par2.zet==0
        par2.t = t;
        par2.k_grid = gr0.k_grid;
    end
    % Interest rate
    par2.R = R_t(t);
    if t==1
        par2.R_last = 1;
    else
        par2.R_last = R_t(t-1);
    end
    
    %%% recreate grid
    gr_now = mainfuncMP_fixed_C(par2,N,[],[],2);

    %%% averages computed using starting distribution!
    % implied total consumption given guess of P (fixed here)
    old_C = 1/(par2.P_old*par2.chi);
    % inner integral for domestic goods
    inner_Cd = gr_now.y_i(:).^((par2.epsi-1)/par2.epsi)'*g0(:);
    % Domestic real output
    Cd = inner_Cd^(par2.epsi/(par2.epsi-1));
    % invert out implied dM to make prices and consumption consistent
    inner_Cf = old_C^((par2.epsi-1)/par2.epsi) - inner_Cd;
    dM_new = inner_Cf / ((Pf.^(-par2.epsi)*(par2.P^par2.epsi)*old_C)^((par2.epsi-1)/par2.epsi));

   
    %%% Measure of firms
    g0  = squeeze(g0);                                 % distribution of active firms
    g0_entrants = pol_t{t}.pr_fstar_enter(:).*gr_now.ps_dist(:);     % distribution of entrants
    g0_stay = pol_t{t}.pr_stay(:).*g0(:);                        % distribution stayers
    g0_exit = (1-pol_t{t}.pr_stay(:)).*g0(:);                    % distribution exiters

    %%% Grids
    incx_gr = pol_t{t}.kpr - (1-par2.delta)*gr_now.k_xgrid;    % gross investment grid
    incx_gr_true = par2.Q*incx_gr.*(incx_gr>0) + ...
        par2.Q*(1-par2.lam)*incx_gr.*(incx_gr<0);          % gross investment grid, factoring in q and Q
    ik_i = incx_gr(:)./gr_now.k_xgrid(:);                   % firm level i/k
    kpr_entrants = pol_t{t}.kpr_up;                          % investment grid of entrants

    %%% Aggregates
    C = 1/(par2.P_old*par2.chi);    % total consumption of final good
    try 
        Cf = cf_i*dM;               % total physical consumption, foreign goods
        P_Cf = Pf*cf_i*dM;          % total real revenue of foreign firms
    catch
        Cf = 0;
        P_Cf = 0;
    end
    P_Cd = gr_now.py_i(:)'*g0(:);   % total real revenue of domestic firms
    Kd = gr_now.k_xgrid(:)'*g0(:);  % total capital stock of domestic firms
    Ld = gr_now.lstar(:)'*g0(:);    % total labor of domestic firms
    TFPR = P_Cd/(Kd^par2.alph * Ld^(1-par2.alph)); % TFPR of domestic firms
    M_norm_t = sum(g0(:)).^(1/(par_t0.epsi-1));
    TFPQ = Cd/(M_norm_t*(Kd^par2.alph * Ld^(1-par2.alph))); % TFPQ of domestic firms

    %%% Investment aggregates (need to account for change in measure)
    window = 0.00;
    % aggregates
    IK = ik_i'*g0_stay(:);                                  % Average I/K of stayers
    I_Gross = incx_gr(:)'*g0_stay(:) + ...
            (1-par2.delta)*(-1*gr_now.k_xgrid(:))'*g0_exit + ...
            kpr_entrants(:)'*g0_entrants(:);                % Gross inv incl stayers, exit and entry
    I_Gross_plus = (incx_gr(:).*(incx_gr(:)>window))'*g0_stay(:) + ...
            kpr_entrants(:)'*g0_entrants(:);                % Gross +ve incl stayers and entry;
    I_Gross_neg = (abs(incx_gr(:)).*(incx_gr(:)<-window))'*g0_stay(:) ...
            + (1-par2.delta)*gr_now.k_xgrid(:)'*g0_exit;      % Gross -ve incl stayers and exits
    I_net = incx_gr_true(:)'*g0_stay(:) + ...
            (1-par2.zet)*par2.Q*(1-par2.lam)*(1-par2.delta)*(-1*gr_now.k_xgrid(:)'*g0_exit) + ...
            par2.Q*kpr_entrants(:)'*g0_entrants(:);                % Net inv incl stayers, exit, and entry
    % Expected values (measures normalized to 1)
    I_Gross_intensive = incx_gr(:)'*g0_stay(:)/sum(g0_stay(:)); % Gross inv along intensive margin
    I_Gross_plus_intensive = (incx_gr(:).*(incx_gr(:)>window))'*g0_stay(:)/sum(g0_stay); % Gross inv along intensive margin
    I_Gross_neg_intensive = (abs(incx_gr(:)).*(incx_gr(:)<-window))'*g0_stay(:)/sum(g0_stay); % Gross inv along intensive margin
    I_Gross_frac_intensive_neg = (incx_gr(:) < -window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of negative investment
    I_Gross_frac_0_intensive = (incx_gr(:)>=-window & incx_gr(:)<=window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of inaction of stayers
    I_Gross_frac_pos_intensive = (incx_gr(:) > window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of pos investment of stayers
    I_Gross_frac_neg_intensive = (incx_gr(:) < -window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of neg investment of stayers

    %%% Other costs, labor and welfare
    operating_cost = pol_t{t}.Ef_pr(:)'*g0_stay(:); % total operating costs paid
    entry_cost = pol_t{t}.Ef_enter(:)'*g0_entrants(:); % total entry costs paid
    L = Ld + I_net + operating_cost + entry_cost + P_Cf; % total labor (domestic + investment goods + fixed costs + export)

    %%% Raw average of TFPQ
    avg_raw_s = gr_now.s_xgrid0(:)'*g0(:)/sum(g0(:));
    
    %%% Raw average of log TFPQ
    avg_log_s = log(gr_now.s_xgrid0(:))'*g0(:)/sum(g0(:));
    sd_log_s = sqrt((log(gr_now.s_xgrid0(:)')-avg_log_s).^2 * g0(:)/sum(g0(:)));
    
    %%% weighted average of TFPQ (PY)
    PY_domestic = gr_now.py_i(:)'*g0(:); % PY as denom
    wgt_py = (gr_now.py_i(:).*g0(:))/PY_domestic; % wgts, PY
    avg_wgt_s = gr_now.s_xgrid0(:)'*wgt_py;
    
    %%% weighted average of TFPQ (Y)
    Y_domestic = gr_now.y_i(:)'*g0(:); % Y as denom
    wgt_y = (gr_now.y_i(:).*g0(:))/Y_domestic; % wgts, Y
    avg_wgt_y_s = gr_now.s_xgrid0(:)'*wgt_y;    
    
    %%% weighted average of TFPQ (K)
    K_domestic = gr_now.k_xgrid(:)'*g0(:); % K as denom
    wgt_k = (gr_now.k_xgrid(:).*g0(:))/K_domestic; % wgts, K
    avg_wgt_k_s = gr_now.s_xgrid0(:)'*wgt_k;    
    
    %%% Average TFPQ of entrants
    avg_raw_s_entrants = gr_now.s_grid(:)'*g0_entrants(:)/sum(g0_entrants(:));
    
    %%% Average TFPQ of exiters
    avg_raw_s_exiters = gr_now.s_xgrid0(:)'*g0_exit(:)/sum(g0_exit(:));
    
    %%% Average TFPQ of continuers
    avg_raw_s_cont = gr_now.s_xgrid0(:)'*g0_stay(:)/sum(g0_stay(:));
    
    %%% MPK dispersion
    mpk_grid = log(gr_now.py_i) - log(gr_now.k_xgrid) + log(par2.alpha_k);
    mean_mpk = mpk_grid(:)'*g0(:)/sum(g0(:));
    sd_mpk = sqrt((mpk_grid(:)'-mean_mpk).^2 * g0(:)/sum(g0(:)));
    
    g0_s = sum(reshape(g0,N.s,N.k),2);
    cdf_g0_s = cumsum(g0_s)/sum(g0_s);
    is = find(cdf_g0_s>=0.5,1,'first');
    s_med = gr_now.s_grid(1:is);
    cdf_med = cdf_g0_s(1:is);
    s_med_2 = gr_now.s_grid(is+1:end);
    cdf_med_2 = cdf_g0_s(is+1:end);
    avg_s_med(t) = s_med(:)'*cdf_med(:)/sum(cdf_med(:));
    avg_s_med_2(t) = s_med_2(:)'*cdf_med_2(:)/sum(cdf_med_2(:));


    %%% store statistics
    C_t(t) = C;
    Cd_t(t) = Cd;
    Cf_t(t) = Cf;
    P_Cd_t(t) = P_Cd;
    P_Cf_t(t) = P_Cf;
    Kd_t(t) = Kd;
    Ld_t(t) = Ld;
    L_t(t) = L;
    TFPR_t(t) = TFPR;
    TFPQ_t(t) = TFPQ;
    Mentrants_t(t) = sum(g0_entrants(:));
    Mactive_t(t) = sum(g0(:));
    Mexit_t(t) = sum(g0_exit(:));
    M_foreign(t) = dM_new;

    IK_t(t) = IK;
    I_Gross_t(t) = I_Gross;
    I_Gross_plus_t(t) = I_Gross_plus;
    I_Gross_neg_t(t) = I_Gross_neg;
    I_net_t(t) = I_net;
    I_Gross_intensive_t(t) = I_Gross_intensive;
    I_Gross_plus_intensive_t(t) = I_Gross_plus_intensive;
	I_Gross_neg_intensive_t(t) = I_Gross_neg_intensive;
    I_Gross_frac_0_intensive_t(t) = I_Gross_frac_0_intensive;
    I_Gross_frac_pos_intensive_t(t) = I_Gross_frac_pos_intensive;
    I_Gross_frac_neg_intensive_t(t) = I_Gross_frac_neg_intensive;
    I_Gross_frac_intensive_neg_t(t) = I_Gross_frac_intensive_neg;

    operating_cost_t(t) = operating_cost;
    entry_cost_t(t) = entry_cost;
    
    avg_raw_s_t(t) = avg_raw_s;
    avg_wgt_s_t(t) = avg_wgt_s;
    avg_wgt_y_s_t(t) = avg_wgt_y_s;
    avg_wgt_k_s_t(t) = avg_wgt_k_s;
    
    avg_raw_s_entrants_t(t) = avg_raw_s_entrants;
    avg_raw_s_exiters_t(t) = avg_raw_s_exiters;
    avg_raw_s_cont_t(t) = avg_raw_s_cont;
    
    mean_mpk_t(t) = mean_mpk;
    sd_mpk_t(t) = sd_mpk;
    sd_log_s_t(t) = sd_log_s;
    
    % average k of exiters
    avg_k_exit_t(t) = gr_now.k_xgrid(:)'*g0_exit(:)/sum(g0_exit(:));
    % average k of continuers
    avg_k_cont_t(t) = gr_now.k_xgrid(:)'*g0_stay(:)/sum(g0_stay(:));
    % average k of entrants
    avg_k_enter_t(t) = kpr_entrants(:)'*g0_entrants(:)/sum(g0_entrants(:));
    % probs k'/k<0 for continuers
    neg_kpr_gr = pol.kpr < gr_now.k_xgrid;
    avg_frac_neg_kpr_t(t) = neg_kpr_gr(:)'*g0_stay(:)/sum(g0_stay(:));

    % push distribution forward
    initval.g0 = g0(:);
    initval.trmat = [];
    output2 = mainfuncMP_fixed_C(par2,N,gr_now,pol_t{t},6,initval);
    g0  = output2.g1(:);

end
x=1;

% back out import penetration
% total_exp = P_t(1:T-1)*C_autarky;
C0_t = 1./(P_t_input(1:T-1)*par2.chi);
total_exp = P_t(1:T-1).*C0_t;
import_pen_t = 1 - P_Cd_t(:)./total_exp(:);

%% compile and pull out stats

out_stats.pol_t = pol_t;
out_stats.gr_t = gr_t;

out_stats.import_pen_t = import_pen_t;
out_stats.sd_mpk_t = sd_mpk_t;
out_stats.TFPQ_t = TFPQ_t;
out_stats.inact_t = I_Gross_frac_0_intensive_t;
out_stats.ipos = I_Gross_frac_pos_intensive_t;
out_stats.ineg = I_Gross_frac_neg_intensive_t;
out_stats.avg_ik = IK_t;
out_stats.agg_ik = I_Gross_t./Kd_t;
out_stats.P_t = P_t(1:T-1);
out_stats.exit_rate_t = Mexit_t(:)./Mactive_t(:);
out_stats.avg_raw_s_t = avg_raw_s_t;

out_stats.avg_wgt_s_t = avg_wgt_s_t;
out_stats.avg_wgt_y_s_t = avg_wgt_y_s_t;
out_stats.avg_wgt_k_s_t = avg_wgt_k_s_t;
out_stats.sd_log_s_t = sd_log_s_t;

%% simulation

 % simulation parameters
    t_burn = 35;
    T_sim = t_burn + T-1;
    Mactive=sum(g0(:));

    %%% Grid policy functions (in index form)
    % continuation kprime policy
    FK_cont = griddedInterpolant({1:N.s,1:N.k},pol_0.i_kpr,'nearest');
    % exit probs, conditioned on state (s,k)
    Fexit = griddedInterpolant({1:N.s,1:N.k},pol_0.pr_stay,'nearest');
    % entry probs, conditioned on state (s)
    Fentry = griddedInterpolant({1:N.s},pol_0.pr_fstar_enter,'nearest');
    % entry kprime policy
    FK_entry = griddedInterpolant({1:N.s},pol_0.i_kpr_up,'nearest');
    % exit choice, conditioned on state (s,k)
    Fstar_exit = griddedInterpolant({1:N.s,1:N.k},pol_0.fstar,'nearest');
    % entry choice, conditioned on state (s)
    Fstar_enter = griddedInterpolant({1:N.s},pol_0.f_star_enter,'nearest');
    % investment policy in levels
    f_kpr = pol_0.kpr;
    f_kpr_up = pol_0.kpr_up;
    
    % actual grid numbers
    [i_s_xgrid , i_k_xgrid] = ndgrid(1:N.s,1:N.k);
    i_s_xgrid = i_s_xgrid(:);
    i_k_xgrid = i_k_xgrid(:);
    
    %%% construct initial panel of active firms
    % frequency at each state (s,fstay,k)
    prk = round(g0_autarky(:)/sum(g0_autarky(:))*10000);
    c_prk = cumsum(prk);
    N_active_star = sum(prk);
    % total length of panel
    N_active = sum(prk); 
    % store initial conditions
    i_k_active = zeros(N_active,1);
    i_s_active = zeros(N_active,1);
    for ii=1:length(prk)
        if ii==1
            istart = 1;
        else
            istart = c_prk(ii-1)+1;
        end
        iend = c_prk(ii);
        
        i_k_active(istart:iend) = repmat(i_k_xgrid(ii),prk(ii),1);
        i_s_active(istart:iend) = repmat(i_s_xgrid(ii),prk(ii),1);
    end
    % invariant distribution of s
    ps_pdf = gr.Ps^1000;
    ps_pdf = ps_pdf(1,:);

    
    i_kcont_sim = cell(T_sim,1);
    i_kstate_stay_sim = cell(T_sim,1);
    i_kenter_sim = cell(T_sim,1);
    i_scont_sim = cell(T_sim,1);
    i_sstate_stay_sim = cell(T_sim,1);
    i_senter_sim = cell(T_sim,1);
    N_active_sim = cell(T_sim,1);
    i_enter_sim = cell(T_sim,1);
    i_stay_sim = cell(T_sim,1);
    i_k_active_sim = cell(T_sim,1);
    i_s_active_sim = cell(T_sim,1);
    i_track_sim = cell(T_sim,1);
    i_s_active_enter_sim = cell(T_sim,1);
    id_sim = cell(T_sim-t_burn+1,1);
    id_allfirms_sim = cell(T_sim-t_burn+1,1);
    year_sim = cell(T_sim-t_burn+1,1);
    year_allfirms_sim = cell(T_sim-t_burn+1,1);
    i_log_mpk_q_sim = cell(T_sim,1);
    i_log_mpk_all_q_sim = cell(T_sim,1);
    log_mpk_i_sim = cell(T_sim,1);
        % track in levels
    kcont_sim = cell(T_sim,1);
    kstate_stay_sim = cell(T_sim,1);
    kenter_sim = cell(T_sim,1);
    k_active_sim = cell(T_sim,1);
    
    % save all initial conditions
    i_kcont_sim{1} = i_k_active;
    i_kenter_sim{1} = [];
    i_scont_sim{1} = i_s_active;
    i_senter_sim{1} = [];
    i_k_active_sim{1} = [];
    i_s_active_sim{1} = [];
    i_s_active_enter_sim{1} = [];
    N_active_sim{1} = N_active;
    kcont_sim{1} = gr.k_grid(i_k_active);

    i_enter_sim{1} = [];
    i_stay_sim{1} = [];
    
    %%% Version 0: Original g(f,s)
    % G(f,s) function
    G_f_s = @(fm,sgrid) fm.*(par_t0.s_sh - (1+sgrid).^par_t0.eta(1) + (1+sgrid).^(par_t0.eta(1)+par_t0.eta(2)));
    % Inverse cdf for sampling f
    igc = @(pr,cm) pr.^(1/par_t0.npow).*cm;
    %%%
    
    tnew = 1;
    sw_track = 0;
    for t=2:T_sim
        
        if t>=t_burn
            %%% Replace policies with transition policies
            % continuation kprime policy
            FK_cont = griddedInterpolant({1:N.s,1:N.k},pol_t{tnew}.i_kpr,'nearest');
            % exit probs, conditioned on state (s,k)
            Fexit = griddedInterpolant({1:N.s,1:N.k},pol_t{tnew}.pr_stay,'nearest');
            % entry probs, conditioned on state (s)
            Fentry = griddedInterpolant({1:N.s},pol_t{tnew}.pr_fstar_enter,'nearest');
            % entry kprime policy
            FK_entry = griddedInterpolant({1:N.s},pol_t{tnew}.i_kpr_up,'nearest');
            % exit choice, conditioned on state (s,k)
            Fstar_exit = griddedInterpolant({1:N.s,1:N.k},pol_t{tnew}.fstar,'nearest');
            % entry choice, conditioned on state (s)
            Fstar_enter = griddedInterpolant({1:N.s},pol_t{tnew}.f_star_enter,'nearest');
            % investment policy in levels
            f_kpr = pol_t{tnew}.kpr;
            f_kpr_up = pol_t{tnew}.kpr_up;
            
            tnew = tnew+1;
            
            % start tracking
            sw_track = 1;
            
            if t==t_burn
                id = 1:N_active_sim{t-1};
                N_track_active = N_active_sim{t-1};
            else
                N_track_active = N_track_active + N_active_sim{t-1};
            end
            
        end

        %%%%%%%%% Current cohort's choices %%%%%%%%%

        N_active = N_active_sim{t-1};
        i_s_active = [i_scont_sim{t-1}(:) ; i_senter_sim{t-1}(:)];
        i_k_active = [i_kcont_sim{t-1}(:) ; i_kenter_sim{t-1}(:)];
        i_track = [0*i_kcont_sim{t-1}(:) ; 1+0*i_kenter_sim{t-1}(:)];
            % track in levels
        k_active = [kcont_sim{t-1}(:) ; kenter_sim{t-1}(:)];

        % draw current period fstay 
        xm = G_f_s(par_t0.mu_fstay_b,gr.s_grid(i_s_active).^((par_t0.epsi-1)/par_t0.epsi));
        f_pot_stay = igc(rand(N_active,1),xm);

        % compute fstar threshold given (s,k)
        f_star_stay = Fstar_exit(i_s_active,i_k_active);

        % make staying choice
        i_stay = f_star_stay > f_pot_stay;

        % keep staying cohort
        Nstay = sum(i_stay);
        i_s_cont = i_s_active(i_stay);
        i_k_cont = i_k_active(i_stay);
            % track in levels
        k_cont = k_active(i_stay);

        % draw next period s for continuers
        i_sprime_stay = sim_discreteMC(gr.Ps,1,0,Nstay,i_s_cont);
        i_sprime_stay = i_sprime_stay(1,:)';

        % compute next period k for continuers
        i_kprime_cont = FK_cont(i_s_cont,i_k_cont);
            % track in levels
        kprime_cont = f_kpr(i_kprime_cont);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%% Entering cohort's choices %%%%%%%%%
        % draw current period fenter
        N_pot_entrants = round((par_t0.Me/Mactive)*N_active_star);

        % draw current period s for potential entrants (iid!)
        i_s_pot_enter = sim_discrete_iid(ps_pdf,N_pot_entrants,'pdf');
        i_s_pot_enter = i_s_pot_enter(:);
        
        % draw current period fstay 
        xm = G_f_s(par_t0.mu_fenter_b,gr.s_grid(i_s_pot_enter).^((par_t0.epsi-1)/par_t0.epsi));
        f_pot_enter = igc(rand(N_pot_entrants,1),xm);

        % compute fstar threshold given this s
        f_star_enter = Fstar_enter(i_s_pot_enter);

        % make entry choice
        i_enter = f_star_enter > f_pot_enter;

        % keep entering cohort
        Nenter = sum(i_enter);
        i_s_enter = i_s_pot_enter(i_enter);

        % draw next period s for entering cohort
        i_sprime_enter = sim_discreteMC(gr.Ps,1,0,Nenter,i_s_enter);
        i_sprime_enter = i_sprime_enter(2,:)';

        % compute next period k for entering cohort
        i_kprime_enter = FK_entry(i_s_enter);
            % track in levels
        kprime_enter = f_kpr_up(i_kprime_enter);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % generate id for continuing and new cohort
        if sw_track
            
            % track all active firms
            id_allfirms_sim{t-t_burn+1} = id(:);
            year_allfirms_sim{t-t_burn+1} = (t-t_burn+1)*ones(N_active,1);
            
            % track continuing firms
            id_stayers = id(i_stay);
            id_sim{t-t_burn+1} = id_stayers(:);
            year_sim{t-t_burn+1} = (t-t_burn+1)*ones(numel(id_stayers),1);
            
            % generate new id for new firms
            id_new = (1:1:Nenter) + N_track_active;
            id = [id_stayers(:) ; id_new(:)];
        end
        
        % compute MPK quantiles using full cross-section
        if t>=t_burn && tnew<17
            
            % pull out y for all firms
            i_y_i = sub2ind([N.s N.k],i_s_active(:),i_k_active(:));
            log_y_i = log(gr_t{tnew-1}.y_i(i_y_i));
            
            % compute MPK and quantiles
            log_mpk_i = log_y_i - (1-par_t0.alph)*log(gr.k_grid(i_k_active(:)));
            n_log_mpk = [-Inf quantile(log_mpk_i,2)+sqrt(eps) Inf];
            
            % pull out quantile indicators for all
            i_log_mpk_all_q = zeros(size(log_mpk_i(:)));
            for iq=1:3
                i_keep = log_mpk_i>n_log_mpk(iq) & log_mpk_i<n_log_mpk(iq+1);
                i_log_mpk_all_q(i_keep) = iq;
            end
            
            % pull out quantile indicators for continuers
            i_log_mpk_q = zeros(size(log_mpk_i(i_stay(:))));
            for iq=1:3
                i_keep = log_mpk_i(i_stay)>n_log_mpk(iq) & log_mpk_i(i_stay)<n_log_mpk(iq+1);
                i_log_mpk_q(i_keep) = iq;
            end
            log_mpk_i_stay = log_mpk_i(i_stay);
    
        else
            i_log_mpk_all_q = nan(size(i_s_active));
            i_log_mpk_q = nan(size(i_s_cont));
            log_mpk_i_stay = nan(size(i_s_cont));
        end
        
        %%% save all variables
        
        % choices 
        i_kcont_sim{t} = i_kprime_cont;
        i_kenter_sim{t} = i_kprime_enter;
        i_scont_sim{t} = i_sprime_stay;
        i_senter_sim{t} = i_sprime_enter;
        N_active_sim{t} = Nstay+Nenter;
        i_enter_sim{t} = i_enter;
        i_stay_sim{t} = i_stay;
        i_kstate_stay_sim{t} = i_k_cont;
        i_sstate_stay_sim{t} = i_s_cont;
            % track in levels
        kcont_sim{t} = kprime_cont;
        kenter_sim{t} = kprime_enter;
        kstate_stay_sim{t} = k_cont;
        
        % states
        i_k_active_sim{t} = i_k_active;
        i_s_active_sim{t} = i_s_active;
        i_track_sim{t} = i_track;
        i_s_active_enter_sim{t} = [i_s_active(i_stay) ; i_s_enter];
            % track in levels
        k_active_sim{t} = k_active;
        
        % mpk
        i_log_mpk_all_q_sim{t} = i_log_mpk_all_q;
        i_log_mpk_q_sim{t} = i_log_mpk_q;
        log_mpk_i_sim{t} = log_mpk_i_stay;
        
        
    end

    % keep only transition data
    i_kcont_sim = i_kcont_sim(t_burn:end-1);
    i_kstate_stay_sim = i_kstate_stay_sim(t_burn:end-1);
    i_scont_sim = i_scont_sim(t_burn:end-1);
    i_sstate_stay_sim = i_sstate_stay_sim(t_burn:end-1);
    i_k_active_sim = i_k_active_sim(t_burn:end-1);
    i_s_active_sim = i_s_active_sim(t_burn:end-1);
    i_s_active_enter_sim = i_s_active_enter_sim(t_burn-1:end-1);
    i_log_mpk_q_sim = i_log_mpk_q_sim(t_burn:end-1);
    i_stay_sim = i_stay_sim(t_burn:end-1);
    i_log_mpk_all_q_sim = i_log_mpk_all_q_sim(t_burn:end-1);
    log_mpk_i_sim = log_mpk_i_sim(t_burn:end-1);
        % track in levels
    kcont_sim = kcont_sim(t_burn:end-1);
    kstate_stay_sim = kstate_stay_sim(t_burn:end-1);
    k_active_sim = k_active_sim(t_burn:end-1);
    
    % inaction window
    inaction_window = .1;
    
    % period 1 data
        % continuers only (investment data)
    out_stats.kpr_vec_t1 = kcont_sim{1};
    out_stats.kstate_vec_t1 = kstate_stay_sim{1}; 
    out_stats.s_vec_t1 = gr.s_grid(i_sstate_stay_sim{1}).^((par_t0.epsi-1)/par_t0.epsi);
    inv_i_t1 = out_stats.kpr_vec_t1 - (1-par_t0.delta)*out_stats.kstate_vec_t1;
    out_stats.inv_i_rate_t1 = inv_i_t1./out_stats.kstate_vec_t1;
    out_stats.inact_i_t1 = abs(out_stats.inv_i_rate_t1)<inaction_window;
    out_stats.i_pos_t1 = out_stats.inv_i_rate_t1>inaction_window;
    out_stats.i_neg_t1 = out_stats.inv_i_rate_t1<-inaction_window;
        % all firms in t=1 (exit)
    out_stats.k_active_sim_t1 = k_active_sim{1};
    out_stats.s_active_sim_t1 = gr.s_grid(i_s_active_sim{1}).^((par_t0.epsi-1)/par_t0.epsi);
    out_stats.i_stay_sim_t1 = i_stay_sim{1};
    
    % convert to single vector
    i_kpr_vec = cell2mat(i_kcont_sim);
    i_kstate_vec = cell2mat(i_kstate_stay_sim);
    i_sstate_vec = cell2mat(i_sstate_stay_sim);
    i_spr_vec = cell2mat(i_scont_sim);
    i_k_vec = cell2mat(i_k_active_sim);
    i_s_vec = cell2mat(i_s_active_sim);
    i_s_lag_vec = cell2mat(i_s_active_enter_sim);
    i_stay_vec = cell2mat(i_stay_sim);
    i_log_mpk_q_sim = cell2mat(i_log_mpk_q_sim);
    i_log_mpk_all_q_sim = cell2mat(i_log_mpk_all_q_sim);
    log_mpk_i_sim = cell2mat(log_mpk_i_sim);
    T2 = T - t_burn;
        % track in levels
    kpr_vec = cell2mat(kcont_sim);
    kstate_vec = cell2mat(kstate_stay_sim);
    k_vec = cell2mat(k_active_sim);
    
    % pull out actual k and s values
    k_active_sim = k_vec;
    s_active_sim = gr.s_grid(i_s_vec).^((par_t0.epsi-1)/par_t0.epsi);
    
    load([file_in transition_file], 'P_Cd_t', 'P_Cf_t')
    import_pen_trend = 1-(P_Cd_t(1:T-1) ./ (P_Cd_t(1:T-1) + P_Cf_t(1:T-1)));

    % ignore panel dimension
    inv_i_t = kpr_vec - (1-par_t0.delta)*kstate_vec;
    inv_i_rate_t = inv_i_t./kstate_vec;
    inact_i_t = abs(inv_i_rate_t)<inaction_window;
    k_i_t = kstate_vec;
    s_i_t = gr.s_grid(i_sstate_vec);
    imp_pen_diff_sample = cell(size(kstate_stay_sim));
    imp_pen_sample = cell(size(kstate_stay_sim));
    for ii=1:length(kstate_stay_sim)
        nsize = numel(kstate_stay_sim{ii});
        imp_pen_diff_sample{ii} = repmat(import_pen_t(ii),nsize,1) - import_pen_trend(ii);
        imp_pen_sample{ii} = repmat(import_pen_t(ii),nsize,1);
        
%         sum(gr.k_grid(i_kstate_stay_sim{ii}))
    end
    
    % only cont
    out_stats.inv_i_t = inv_i_t;
    out_stats.inact_i_t = inact_i_t;
    out_stats.k_i_t = k_i_t;
    out_stats.s_i_t = s_i_t;
    out_stats.imp_pen_diff_sample = cell2mat(imp_pen_diff_sample);
    out_stats.imp_pen_sample = cell2mat(imp_pen_sample);
    out_stats.id_sim = cell2mat(id_sim(1:end-1));
    out_stats.year_sim = cell2mat(year_sim(1:end-1));
    out_stats.i_log_mpk_q_sim = i_log_mpk_q_sim;
    out_stats.log_mpk_i_sim = log_mpk_i_sim;
    
    % incl cont / exit
    out_stats.i_stay_vec = i_stay_vec;
    out_stats.k_active_sim = k_active_sim;
    out_stats.s_active_sim = s_active_sim;
    out_stats.i_log_mpk_all_q_sim = i_log_mpk_all_q_sim;
    clear imp_pen_diff_sample
    imp_pen_diff_sample = cell(size(i_stay_sim));
    imp_pen_sample = cell(size(i_stay_sim));
    for ii=1:length(i_stay_sim)
        nsize = numel(i_stay_sim{ii});
        imp_pen_diff_sample{ii} = repmat(import_pen_t(ii),nsize,1) - import_pen_trend(ii);
        imp_pen_sample{ii} = repmat(import_pen_t(ii),nsize,1);
    end
    out_stats.imp_pen_diff_sample_2 = cell2mat(imp_pen_diff_sample);
    out_stats.imp_pen_sample_2 = cell2mat(imp_pen_sample);
    out_stats.year_allfirms_sim = cell2mat(year_allfirms_sim(1:end-1));
    out_stats.id_allfirms_sim = cell2mat(id_allfirms_sim(1:end-1));
    
%     for ii=1:length(i_kstate_stay_sim)
%         inv_i_t = gr.k_grid(i_kcont_sim) - (1-par_t0.delta)*gr.k_grid(i_kstate_stay_sim);
%         inact_i_t = abs(inv_i_t)<.01;
%     end

