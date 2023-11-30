function p_stats = g_event_study(init_ss,trans_x)

T = trans_x.T;

% load initial steady state
load(init_ss, 'Px', 'par2', 'gr', 'N', 'Pnow', 'old_C','pol')
g0_autarky = squeeze(Px.g_ss);
P_autarky = Pnow;
C_autarky = old_C;
par_t0 = par2;
kmax_autarky = gr.k_grid(end);
par_t0.kstar = kmax_autarky;
par_t0.P0 = P_autarky;
par_t0.C0 = C_autarky;
pol_0 = pol;

%% Run simulation

% simulation parameters
t_burn = 50;
T_sim = t_burn + T-1;
Mactive=sum(g0_autarky(:));

% price / agg C at autarky
Px = P_autarky;
Yd = C_autarky;

%%% Grid policy functions (in linear form)
% continuation kprime policy
FK_cont = griddedInterpolant({1:N.s,gr.k_grid},pol_0.kpr,'linear');
% exit probs, conditioned on state (s,k)
Fexit = griddedInterpolant({1:N.s,gr.k_grid},pol_0.pr_stay,'linear');
% entry probs, conditioned on state (s)
Fentry = griddedInterpolant({1:N.s},pol_0.pr_fstar_enter,'linear');
% entry kprime policy
FK_entry = griddedInterpolant({1:N.s},pol_0.kpr_up,'linear');
% exit choice, conditioned on state (s,k)
Fstar_exit = griddedInterpolant({1:N.s,gr.k_grid},pol_0.fstar,'linear');
% entry choice, conditioned on state (s)
Fstar_enter = griddedInterpolant({1:N.s},pol_0.f_star_enter,'linear');

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


kcont_sim = cell(T_sim,1);
kstate_stay_sim = cell(T_sim,1);
kenter_sim = cell(T_sim,1);
i_scont_sim = cell(T_sim,1);
i_sstate_stay_sim = cell(T_sim,1);
i_senter_sim = cell(T_sim,1);
N_active_sim = cell(T_sim,1);
i_enter_sim = cell(T_sim,1);
i_stay_sim = cell(T_sim,1);
k_active_sim = cell(T_sim,1);
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

% save all initial conditions
kcont_sim{1} = gr.k_grid(i_k_active);
kenter_sim{1} = [];
i_scont_sim{1} = i_s_active;
i_senter_sim{1} = [];
k_active_sim{1} = [];
i_s_active_sim{1} = [];
i_s_active_enter_sim{1} = [];
N_active_sim{1} = N_active;

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
    k_active = [kcont_sim{t-1}(:) ; kenter_sim{t-1}(:)];
    i_track = [0*kcont_sim{t-1}(:) ; 1+0*kenter_sim{t-1}(:)];

    % draw current period fstay 
    xm = G_f_s(par_t0.mu_fstay_b,gr.s_grid(i_s_active).^((par_t0.epsi-1)/par_t0.epsi));
    f_pot_stay = igc(rand(N_active,1),xm);

    % compute fstar threshold given (s,k)
    f_star_stay = Fstar_exit(i_s_active,k_active);

    % make staying choice
    i_stay = f_star_stay > f_pot_stay;

    % keep staying cohort
    Nstay = sum(i_stay);
    i_s_cont = i_s_active(i_stay);
    k_cont = k_active(i_stay);

    % draw next period s for continuers
    i_sprime_stay = sim_discreteMC(gr.Ps,1,0,Nstay,i_s_cont);
    i_sprime_stay = i_sprime_stay(1,:)';

    % compute next period k for continuers
    kprime_cont = FK_cont(i_s_cont,k_cont);
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
    kprime_enter = FK_entry(i_s_enter);

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

    %%%% need to recode because old version uses on the grid

    % compute MPK quantiles using full cross-section
    if t>=t_burn && tnew<T+1

        % labor given s, k, P, C
        Afac = par_t0.A*Px*Yd^(1/par_t0.epsi);
        s_grid_rev = gr.s_grid.^((par_t0.epsi-1)/par_t0.epsi);
        lstar = (par_t0.alpha_l*Afac.*s_grid_rev(i_s_active).*k_active.^par_t0.alpha_k).^(1/(1-par_t0.alpha_l));    

        % revenue
        y_i = Afac.*s_grid_rev(i_s_active).*k_active.^par_t0.alpha_k.*lstar.^par_t0.alpha_l;
        log_y_i = log(y_i);

        % compute MPK and quantiles
        log_mpk_i = log_y_i - (1-par_t0.alph)*log(k_active(:));
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
        
        % pull out mpk for continuers
        log_mpk_i_temp = log_mpk_i(i_stay(:));
        clear log_mpk_i
        log_mpk_i = log_mpk_i_temp;
        clear log_mpk_i_temp

    else
        i_log_mpk_all_q = nan(size(i_s_active));
        i_log_mpk_q = nan(size(i_s_cont));
        log_mpk_i = nan(size(i_s_cont));
    end
    %%%%%%%%%%%%%%%%

    %%% save all variables

    % choices 
    kcont_sim{t} = kprime_cont;
    kenter_sim{t} = kprime_enter;
    i_scont_sim{t} = i_sprime_stay;
    i_senter_sim{t} = i_sprime_enter;
    N_active_sim{t} = Nstay+Nenter;
    i_enter_sim{t} = i_enter;
    i_stay_sim{t} = i_stay;
    kstate_stay_sim{t} = k_cont;
    i_sstate_stay_sim{t} = i_s_cont;

    % states
    k_active_sim{t} = k_active;
    i_s_active_sim{t} = i_s_active;
    i_track_sim{t} = i_track;
    i_s_active_enter_sim{t} = [i_s_active(i_stay) ; i_s_enter];

    % mpk
    i_log_mpk_all_q_sim{t} = i_log_mpk_all_q;
    i_log_mpk_q_sim{t} = i_log_mpk_q;
    log_mpk_i_sim{t} = log_mpk_i;


end

% keep selected data
kcont_sim = kcont_sim(t_burn:end-1);
kstate_stay_sim = kstate_stay_sim(t_burn:end-1);
i_sstate_stay_sim = i_sstate_stay_sim(t_burn:end-1);
i_log_mpk_q_sim = i_log_mpk_q_sim(t_burn:end-1);
log_mpk_i_sim = log_mpk_i_sim(t_burn:end-1);

% convert to single vector
kpr_vec = cell2mat(kcont_sim);
kstate_vec = cell2mat(kstate_stay_sim);
i_sstate_vec = cell2mat(i_sstate_stay_sim);
i_log_mpk_q_sim = cell2mat(i_log_mpk_q_sim);
log_mpk_i_sim = cell2mat(log_mpk_i_sim);

% ignore panel dimension
inv_i_t = kpr_vec - (1-par_t0.delta)*kstate_vec;
inv_i_rate_t = inv_i_t./kstate_vec;
inact_i_t = abs(inv_i_rate_t)<.1;
k_i_t = kstate_vec;
s_i_t = gr.s_grid(i_sstate_vec);

% export as a structure
p_stats.inv_i_t = inv_i_t;
p_stats.inact_i_t = inact_i_t;
p_stats.k_i_t = k_i_t;
p_stats.s_i_t = s_i_t;
p_stats.id_sim = cell2mat(id_sim(1:end-1));
p_stats.year_sim = cell2mat(year_sim(1:end-1));
p_stats.i_log_mpk_q_sim = i_log_mpk_q_sim;
p_stats.log_mpk_i_sim = log_mpk_i_sim;


    

