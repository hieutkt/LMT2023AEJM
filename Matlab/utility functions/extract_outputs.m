function extract_outputs(dat_path,file0,file1,file2,file0_q1,file1_q1,file2_q1,versionx)
% extract_outputs.m
%
%   Use this function to load the results from the raw workspaces, and
%   compact into a single structure that saves the desired outputs into a
%   new workspace called ws_plots

clc

% store current path
pathnow = pwd;

cd(dat_path)

%% q<Q

%%% autarky
load(file0,'output_moment', 'par2', 'pol', 'gr', 'Px', 'N')
stats_autarky_qb = output_moment.stats;
avg_autarky = output_moment.avg;
moments_qb = output_moment.moments;
bhat_probit_qb = output_moment.bhat_probit;
avg_autarky_qb = avg_autarky;
par_t0 = par2;
par_qb = par2;
pol_qb = pol;
gr_qb = gr;
Px_qb = Px;
N_qb = N;
output_moment_qb_t0 = output_moment;
clear output_moment par2 Px gr pol

%%% trade
load(file1, 'output_moment', 'Pf', 'par2', 'gr')
stats_trade_qb = output_moment.stats;
avg_trade = output_moment.avg;
avg_trade_qb = output_moment.avg;
Pf_qb = Pf;
par_q1 = par2;
gr_qb_trade = gr;
output_moment_qb_trade = output_moment;
clear output_moment Pf par2 gr

%%% transition path
load(file2, '*_t', 'T', 'pol_t')
pol_t_qb = pol_t;
clear pol_t

%%% compute s-s comparisons
Y_manu_qb = avg_trade.Cd/avg_autarky.C; 
C_all_qb = avg_trade.C/avg_autarky.C; 
K_qb = avg_trade.Kd/avg_autarky.Kd; % capital
L_manu_qb = avg_trade.Ld/avg_autarky.Ld; % manufacturing only
tfp_qb = avg_trade.TFPR/avg_autarky.TFPR; % TFPR
L_all_qb = avg_trade.L/avg_autarky.L; % total labor
P_qb = avg_trade.P/avg_autarky.P;
 
% pack into a vector
welfare_ss_qb = exp((avg_trade.Uss-avg_autarky.Uss)*(1-par_t0.betaa))-1;
stats_qb = [Y_manu_qb C_all_qb K_qb L_manu_qb tfp_qb L_all_qb P_qb];

% save as "qb" vars
P_plot_qb = [avg_autarky.P P_t(1:end-1)'];
tfp_t_qb = [avg_autarky.TFPR TFPR_t(:)'];
l_t_qb = [avg_autarky.Ld Ld_t(:)'] ;
l_all_t_qb = [avg_autarky.L L_t(:)'];
inv_net_t_qb = [avg_autarky.I_net I_net_t(:)'];
k_t_qb = [avg_autarky.Kd Kd_t(:)'];
kstay_t_qb = [avg_autarky.Kd_stay Kd_t(:)'];
C_t_qb = [avg_autarky.C C_t(:)'];
Cd_t_qb = [avg_autarky.C Cd_t(:)'];
Cf_t_qb = [avg_autarky.Cf Cf_t(:)'];
P_Cd_t_qb = [avg_autarky.P_Cd P_Cd_t(:)'];
P_Cf_t_qb = [avg_autarky.P_Cf P_Cf_t(:)'];
tfpq_t_qb = tfp_t_qb(:)./P_plot_qb(:);
IK_t_qb = [avg_autarky.IK IK_t(:)'];
% oc_t_qb = oc_t;
% ec_t_qb = ec_t;
MM_t_qb = [avg_autarky.Mactive Mactive_t(:)'];
Mentry_t_qb = [avg_autarky.Mentrants Mentrants_t(:)'];
Mexit_t_qb = [avg_autarky.Mexit Mexit_t(:)'];

inc_Q_t_qb = [avg_autarky.I_Gross I_Gross_t(:)'];
inc_Q_plus_t_qb = [avg_autarky.I_Gross_plus I_Gross_plus_t(:)'];
inc_Q_neg_t_qb = [avg_autarky.I_Gross_neg I_Gross_neg_t(:)'];
inc_Q_0_t_qb = [avg_autarky.I_Gross_frac_0_intensive I_Gross_frac_0_intensive_t(:)'];
inc_Q_pos_10p_qb = [avg_autarky.I_Gross_frac_pos_intensive I_Gross_frac_pos_intensive_t(:)'];
inc_Q_neg_10p_qb = [avg_autarky.I_Gross_frac_neg_intensive I_Gross_frac_neg_intensive_t(:)'];
inc_Q_0_10p_qb = 1 - inc_Q_pos_10p_qb - inc_Q_neg_10p_qb;

% back out mass of stayers
Mstay = avg_autarky.Mactive - avg_autarky.Mexit;
Mstay_t_qb = [Mstay Mactive_t(:)' - Mexit_t(:)'];

% intensive margin
inv_intensive_t_qb = [avg_autarky.I_Gross_intensive I_Gross_intensive_t(:)'];
inv_neg_intensive_t_qb = [avg_autarky.I_Gross_neg_intensive I_Gross_neg_intensive_t(:)'];
inv_pos_intensive_t_qb = [avg_autarky.I_Gross_plus_intensive I_Gross_plus_intensive_t(:)'];
I_exit_t = abs(I_Gross_neg_t(:)') - abs(I_Gross_neg_intensive_t(:)'.*Mstay_t_qb(2:end));
I_exit = abs(avg_autarky.I_Gross_neg) - abs(avg_autarky.I_Gross_neg_intensive*Mstay);
I_exit_t_qb = [I_exit I_exit_t(:)'];

% TFP
avg_raw_s_t_qb = [avg_autarky.avg_raw_s avg_raw_s_t(:)'];
avg_wgt_s_t_qb = [avg_autarky.avg_wgt_s avg_wgt_s_t(:)'];
avg_wgt_y_s_t_qb = [avg_autarky.avg_wgt_y_s avg_wgt_y_s_t(:)'];
try
    avg_wgt_k_s_t_qb = [avg_autarky.avg_wgt_k_s avg_wgt_k_s_t(:)'];
catch
    disp('Old version, no k weighted TFPR')
end

% MPK
sdmpk_t_qb = [stats_autarky_qb.sd_mpk sd_mpk_t(:)'];

% Physical TFPQ: Domestic physical out / factors*scale
M_norm_t = MM_t_qb.^(1/(par_t0.epsi-1));
TFPQ_t_qb = Cd_t_qb./(M_norm_t.*(k_t_qb.^par_t0.alph .* l_t_qb.^(1-par_t0.alph)));
TFPQ_uncorrected_t_qb = Cd_t_qb./(k_t_qb.^par_t0.alph .* l_t_qb.^(1-par_t0.alph));
M_norm_t_qb = M_norm_t;

% welfare gains
V=0;
diffv = 1e-7;
U_vec_qb = zeros(1000,1);
for ii=1:1000
    
    if ii<=T
        C_now = C_t_qb(ii);
        L_now = l_all_t_qb(ii);
    else
        C_now = avg_trade.C;
        L_now = avg_trade.L;
    end
    U_now = log(C_now) - par_t0.chi*L_now;
    V_now = V + par_t0.betaa^(ii-1)*U_now;
    if abs(V_now-V)<diffv
        break
    else
        V = V_now;
    end
    
    U_vec_qb(ii) = U_now;
    
end
U_shock_qb = V;
welfare_shock_qb = exp((U_shock_qb-avg_autarky.Uss)*(1-par_t0.betaa))-1;
U_vec_qb = U_vec_qb(1:ii-1);

%%% other statistics from steady-state
% capital weighted exit rates
wgt_k = Px_qb.g_ss(:).*gr_qb.k_xgrid(:)/(Px_qb.g_ss(:)'*gr_qb.k_xgrid(:));
avg_wgt_k_exit_qb = 1-wgt_k(:)'*pol_qb.pr_stay(:);

clear avg_autarky avg_trade

%% q=Q

try

    %%% autarky
    load(file0_q1, 'output_moment', 'pol', 'gr', 'Px', 'par2', 'N')
    stats_autarky_q1 = output_moment.stats;
    avg_autarky = output_moment.avg;
    avg_autarky_q1 = avg_autarky;
    output_moment_q1_t0 = output_moment;
    pol_q1 = pol;
    gr_q1 = gr;
    Px_q1 = Px;
    par_q1 = par2;
    N_q1 = N;
    par_t0 = par2;
    clear output_moment pol gr Px par2 N

    %%% trade
    load(file1_q1, 'output_moment', 'Pf', 'dM', 'gr')
    avg_trade = output_moment.avg;
    avg_trade_q1 = output_moment.avg;
    Pf_q1 = Pf;
    dM_q1 = dM;
    gr_q1_trade = gr;
    output_moment_q1_trade = output_moment;
    clear output_moment Pf gr

    %%% transition path
    load(file2_q1, '*_t', 'T', 'pol_t')
    pol_t_q1 = pol_t;
    clear pol_t

    %%% compute s-s comparisons
    Y_manu_q1 = avg_trade.Cd/avg_autarky.C; 
    C_all_q1 = avg_trade.C/avg_autarky.C; 
    K_q1 = avg_trade.Kd/avg_autarky.Kd; % capital
    L_manu_q1 = avg_trade.Ld/avg_autarky.Ld; % manufacturing only
    tfp_q1 = avg_trade.TFPR/avg_autarky.TFPR; % TFPR
    L_all_q1 = avg_trade.L/avg_autarky.L; % total labor
    P_q1 = avg_trade.P/avg_autarky.P;

    % pack into a vector
    welfare_ss_q1 = exp((avg_trade.Uss-avg_autarky.Uss)*(1-par_t0.betaa))-1;
    stats_q1 = [Y_manu_q1 C_all_q1 K_q1 L_manu_q1 tfp_q1 L_all_q1 P_q1];

    % save as "qb" vars
    P_plot_q1 = [avg_autarky.P P_t(1:end-1)'];
    tfp_t_q1 = [avg_autarky.TFPR TFPR_t(:)'];
    l_t_q1 = [avg_autarky.Ld Ld_t(:)'] ;
    l_all_t_q1 = [avg_autarky.L L_t(:)'];
    inv_net_t_q1 = [avg_autarky.I_net I_net_t(:)'];
    k_t_q1 = [avg_autarky.Kd Kd_t(:)'];
    C_t_q1 = [avg_autarky.C C_t(:)'];
    Cd_t_q1 = [avg_autarky.C Cd_t(:)'];
    Cf_t_q1 = [avg_autarky.Cf Cf_t(:)'];
    P_Cd_t_q1 = [avg_autarky.P_Cd P_Cd_t(:)'];
    P_Cf_t_q1 = [avg_autarky.P_Cf P_Cf_t(:)'];
    tfpq_t_q1 = tfp_t_q1(:)./P_plot_q1(:);
    IK_t_q1 = [avg_autarky.IK IK_t(:)'];
    % oc_t_q1 = oc_t;
    % ec_t_q1 = ec_t;
    MM_t_q1 = [avg_autarky.Mactive Mactive_t(:)'];
    Mentry_t_q1 = [avg_autarky.Mentrants Mentrants_t(:)'];
    Mexit_t_q1 = [avg_autarky.Mexit Mexit_t(:)'];

    inc_Q_t_q1 = [avg_autarky.I_Gross I_Gross_t(:)'];
    inc_Q_plus_t_q1 = [avg_autarky.I_Gross_plus I_Gross_plus_t(:)'];
    inc_Q_neg_t_q1 = [avg_autarky.I_Gross_neg I_Gross_neg_t(:)'];
    inc_Q_0_t_q1 = [avg_autarky.I_Gross_frac_0_intensive I_Gross_frac_0_intensive_t(:)'];

    % back out mass of stayers
    Mstay = avg_autarky.Mactive - avg_autarky.Mexit;
    Mstay_t_q1 = [Mstay Mactive_t(:)' - Mexit_t(:)'];

    % intensive margin
    inv_intensive_t_q1 = [avg_autarky.I_Gross_intensive I_Gross_intensive_t(:)'];
    inv_neg_intensive_t_q1 = [avg_autarky.I_Gross_neg_intensive I_Gross_neg_intensive_t(:)'];
    inv_pos_intensive_t_q1 = [avg_autarky.I_Gross_plus_intensive I_Gross_plus_intensive_t(:)'];
    I_exit_t = abs(I_Gross_neg_t(:)') - abs(I_Gross_neg_intensive_t(:)'.*Mstay_t_q1(2:end));
    I_exit = abs(avg_autarky.I_Gross_neg) - abs(avg_autarky.I_Gross_neg_intensive*Mstay);
    I_exit_t_q1 = [I_exit I_exit_t];

    % TFP
    avg_raw_s_t_q1 = [avg_autarky.avg_raw_s avg_raw_s_t(:)'];
    avg_wgt_s_t_q1 = [avg_autarky.avg_wgt_s avg_wgt_s_t(:)'];
    avg_wgt_y_s_t_q1 = [avg_autarky.avg_wgt_y_s avg_wgt_y_s_t(:)'];
    try
        avg_wgt_k_s_t_q1 = [avg_autarky.avg_wgt_k_s avg_wgt_k_s_t(:)'];
    catch
        disp('Old version, no k weighted TFPR')
    end

    % MPK
    sdmpk_t_q1 = [stats_autarky_q1.sd_mpk sd_mpk_t(:)'];
    
    % Physical TFPQ: Domestic physical out / factors*scale
    M_norm_t = MM_t_q1.^(1/(par_t0.epsi-1));
    TFPQ_t_q1 = Cd_t_q1./(M_norm_t.*(k_t_q1.^par_t0.alph .* l_t_q1.^(1-par_t0.alph)));
    TFPQ_t_q1(2) = TFPQ_t_q1(1);
    TFPQ_uncorrected_t_q1 = Cd_t_q1./(k_t_q1.^par_t0.alph .* l_t_q1.^(1-par_t0.alph));
    M_norm_t_q1 = M_norm_t;

    V=0;
    diffv = 1e-7;
    U_vec_q1 = zeros(1000,1);
    for ii=1:1000

        if ii<=T
            C_now = C_t_q1(ii);
            L_now = l_all_t_q1(ii);
        else
            C_now = avg_trade.C;
            L_now = avg_trade.L;
        end
        U_now = log(C_now) - par_t0.chi*L_now;
        V_now = V + par_t0.betaa^(ii-1)*U_now;
        if abs(V_now-V)<diffv
            break
        else
            V = V_now;
        end

        U_vec_q1(ii) = U_now;

    end
    U_shock_q1 = V;
    welfare_shock_q1 = exp((U_shock_q1-avg_autarky.Uss)*(1-par_t0.betaa))-1;
    U_vec_q1 = U_vec_q1(1:ii-1);

    clear avg_autarky avg_trade

catch

    disp('q1 not supplied')

end

%% Compute thresholds and (S,s) bands

try
    % evaluate at mean f
    try
        switch versionx
            case('normal')
                fx=exp(par_qb.mu_fstay + 1/2*par_qb.sigma_ufstay^2);

            case('uniform')
                if par_qb.npow>0
                    ffx = @(y) (1/(par_qb.npow+1))*y.^(par_qb.npow+1);
                elseif par_qb.npow==0
                    ffx = @(y) (y.*log(y) - y);
                end
                fx=(ffx(par_qb.mu_fstay_b) - ffx(par_qb.mu_fstay_a)) / (par_qb.mu_fstay_b-par_qb.mu_fstay_a);
        end

    catch
        disp('SUPPLY VERSION!')
        return

    end
    out_qb = thresholds(par_qb, N_qb, pol_qb, gr_qb, [], output_moment_qb_t0, [], fx);
    out_qb_trade = thresholds(par_qb, N_qb, pol_t_qb{1}, gr_qb, [], output_moment_qb_trade, [], fx);
    out_q1 = thresholds(par_q1, N_q1, pol_q1, gr_q1, [], output_moment_q1_t0, [], fx);
    out_q1_trade = thresholds(par_q1, N_q1, pol_t_q1{1}, gr_q1, [], output_moment_q1_trade, [], fx);

    % close all
    % nsig= 5;
    % out_qb_s = cell(nsig,1);
    % out_q1_s = cell(nsig,1);
    % hold on
    % for ii=1:nsig
    %     % evaluate at mean f
    % fx=exp(par_qb.mu_fstay + 1/2*par_qb.sigma_ufstay^2 + .25*ii*par_qb.sigma_ufstay);
    % out_qb_s{ii} = thresholds(par_qb, N_qb, pol_qb, gr_qb, [], [], [], fx);
    % out_q1_s{ii} = thresholds(par_q1, N_q1, pol_q1, gr_q1, [], [], [], fx);
    % 
    % plot(gr_qb.k_grid,out_qb_s{ii}.s_fx_star,'linewidth',3);
    % end
    % 
    % hold off

catch

    disp('not computing thresholds')

end

%% Export statistics

try

    close all
    clc

    % table 5 main stats
    tab_stats_qb = ([avg_trade_qb.C/avg_autarky_qb.C ...
        avg_trade_qb.Kd/avg_autarky_qb.Kd ...
        avg_trade_qb.Ld/avg_autarky_qb.Ld ...
        avg_trade_qb.Mactive/avg_autarky_qb.Mactive ...
        avg_trade_qb.avg_raw_s/avg_autarky_qb.avg_raw_s]-1)*100;


    % table 6 (old)
    tab_old.labels = {'sd mpk', 'rho_hat'};
    tab_old.qb = [stats_autarky_qb.sd_mpk stats_autarky_qb.rho_hat_condn(2:end)'];
    tab_old.q1 = [stats_autarky_q1.sd_mpk stats_autarky_q1.rho_hat_condn(2:end)'];
    disp('q<Q, table 6')
    disp(tab_old.qb)
    disp('q=Q, table 6')
    disp(tab_old.q1)

    % table 6 (new)
    tab_new.labels = {'sd mpk', 'pr(1,1)', 'pr(3,3)', 'pr(3,3)'};
    tab_new.qb = [stats_autarky_qb.sd_mpk stats_autarky_qb.pr_mpk_raw_norm(1,1) stats_autarky_qb.pr_mpk_raw_norm(2,2) stats_autarky_qb.pr_mpk_raw_norm(3,3)];
    tab_new.q1 = [stats_autarky_q1.sd_mpk stats_autarky_q1.pr_mpk_raw_norm(1,1) stats_autarky_q1.pr_mpk_raw_norm(2,2) stats_autarky_q1.pr_mpk_raw_norm(3,3)];
    disp('Pr(x(t+1)|x(t), q<Q')
    disp(round([[1;2;3] stats_autarky_qb.pr_mpk_raw_norm],2))
    disp('Pr(x(t+1)|x(t), q=Q')
    disp(round([[1;2;3] stats_autarky_q1.pr_mpk_raw_norm],2))

    % table 7
    tab_closed = [avg_autarky_qb.C avg_autarky_qb.Kd avg_autarky_qb.Ld avg_autarky_qb.Mactive avg_autarky_qb.avg_raw_s];
    tab_open = [avg_trade_qb.C avg_trade_qb.Kd avg_trade_qb.Ld avg_trade_qb.Mactive avg_trade_qb.avg_raw_s];
    disp('C, K, L all, L manu, M, TFPQ')
    disp('closed')
    disp(tab_closed)
    disp('open')
    disp(tab_open)
    disp('open/closed, in percent')
    disp(100*(tab_open./tab_closed-1))

    % table welfare
    ws_plots.tab_welfare = [welfare_ss_q1 welfare_shock_q1 ; welfare_ss_qb welfare_shock_qb];

    %%% calibration outcomes 
    try
        % load full set of targets
        load('/Users/EugeneTan/Dropbox/1. Research/1 Active projects/006 Peru/Andrea_Pamela/Model/extensive_margin/_baseline_scripts/model_fixed_P_uniform/cali/targets_5p_1sig', 'dat_in')
        dat_moments = dat_in;
        xx = zeros(size(dat_moments));
        xx(1) = output_moment_qb_t0.moments.I_Gross_frac_intensive_neg;
        xx(2) = output_moment_qb_t0.avg.Mexit/output_moment_qb_t0.avg.Mactive;
        xx(3) = output_moment_qb_t0.moments.avg_k_exit/output_moment_qb_t0.moments.avg_k_all;
        xx(4) = output_moment_qb_t0.avg.avg_raw_s_exiters/output_moment_qb_t0.avg.avg_raw_s;
        xx(5) = output_moment_qb_t0.stats.rho_log_tfpr;
        xx(6) = output_moment_qb_t0.stats.sd_log_tfpr;
        xx(7) = output_moment_qb_t0.bhat_probit(3) / output_moment_qb_t0.bhat_probit(2);
        stats_out.moments_outcome = [xx(:) dat_moments(:)];
    catch
        xx = zeros(8,1);
        xx(1) = output_moment_qb_t0.moments.I_Gross_frac_intensive_neg;
        xx(5) = output_moment_qb_t0.avg.Mexit/output_moment_qb_t0.avg.Mactive;
        xx(6) = output_moment_qb_t0.moments.avg_k_exit/output_moment_qb_t0.moments.avg_k_all;
        xx(7) = output_moment_qb_t0.avg.avg_raw_s_exiters/output_moment_qb_t0.avg.avg_raw_s;
        xx(3) = output_moment_qb_t0.stats.rho_log_tfpr;
        xx(4) = 0.75*output_moment_qb_t0.stats.sd_log_tfpr;
        xx(2) = output_moment_qb_t0.bhat_probit(3) / output_moment_qb_t0.bhat_probit(2);
        xx(8) = output_moment_qb_t0.stats.sd_ik;
        stats_out.moments_outcome = xx(:);
    end


    % untargetted
    stats_out.model.pr_mpk = stats_autarky_qb.pr_mpk_raw_norm;
    stats_out.model.sd_mpk = stats_autarky_qb.sd_mpk;
    stats_out.data.sd_mpk = 1.47;

catch

    disp('q1 model not supplied')

end
%% Construct workspace for time series

try

    % prices
    ws_plots.P_plot_qb = P_plot_qb;
    ws_plots.P_plot_q1 = P_plot_q1;

    % cons
    ws_plots.C_t_q1 = C_t_q1;
    ws_plots.C_t_qb = C_t_qb;
    
    % domestic production
    ws_plots.Cd_t_q1 = Cd_t_q1;
    ws_plots.Cd_t_qb = Cd_t_qb;

    % capital
    ws_plots.k_t_q1 = k_t_q1;
    ws_plots.k_t_qb = k_t_qb;

    % Labor manu
    ws_plots.l_t_q1 = l_t_q1;
    ws_plots.l_t_qb = l_t_qb;
    
    % Labor all
    ws_plots.l_all_t_qb = l_all_t_qb;
    ws_plots.l_all_t_q1 = l_all_t_q1;

    % startup rate
    ws_plots.Mentry_t_qb = Mentry_t_qb;
    ws_plots.Mentry_t_q1 = Mentry_t_q1;

    % exit rate
    ws_plots.Mexit_t_qb = Mexit_t_qb;
    ws_plots.Mexit_t_q1 = Mexit_t_q1;

    % Mass of active firms
    ws_plots.MM_t_q1 = MM_t_q1;
    ws_plots.MM_t_qb = MM_t_qb;

    % Mass of continuing firms
    ws_plots.MM_cont_t_q1 = Mstay_t_q1;
    ws_plots.MM_cont_t_qb = Mstay_t_qb;

    % Investment
    ws_plots.inc_Q_t_qb = inc_Q_t_qb;
    ws_plots.inc_Q_t_q1 = inc_Q_t_q1;

    % Import penetration
    ws_plots.imp_t_qb = P_Cf_t_qb(:)./(P_Cf_t_qb(:)+P_Cd_t_qb(:));
    ws_plots.imp_t_q1 = P_Cf_t_q1(:)./(P_Cf_t_q1(:)+P_Cd_t_q1(:));

    % Imports relative to domestic output
    ws_plots.rel_Cf_Cd_qb = P_Cf_t_qb(:)./P_Cd_t_qb(:);
    ws_plots.rel_Cf_Cd_q1 = P_Cf_t_q1(:)./P_Cd_t_q1(:);

    % TFP (raw average)
    ws_plots.avg_raw_s_t_q1 = avg_raw_s_t_q1;
    ws_plots.avg_raw_s_t_qb = avg_raw_s_t_qb;

    % TFP (PY weighted average)
    ws_plots.avg_wgt_s_t_q1 = avg_wgt_s_t_q1;
    ws_plots.avg_wgt_s_t_qb = avg_wgt_s_t_qb;

    % TFP (Y weighted average)
    ws_plots.avg_wgt_y_s_t_q1 = avg_wgt_y_s_t_q1;
    ws_plots.avg_wgt_y_s_t_qb = avg_wgt_y_s_t_qb;

    % TFP (K weighted average)
    ws_plots.avg_wgt_k_s_t_q1 = avg_wgt_k_s_t_q1;
    ws_plots.avg_wgt_k_s_t_qb = avg_wgt_k_s_t_qb;

    % Measured TFPQ
    ws_plots.tfpq_t_q1 = tfpq_t_q1;
    ws_plots.tfpq_t_qb = tfpq_t_qb;

    % Measured TFPR
    ws_plots.tfpr_t_q1 = tfp_t_q1;
    ws_plots.tfpr_t_qb = tfp_t_qb;
    
    % True physical TFPQ
    ws_plots.TFPQ_t_qb = TFPQ_t_qb;
    ws_plots.TFPQ_t_q1 = TFPQ_t_q1;
    
    % Uncorrected TFOQ
    ws_plots.TFPQ_uncorrected_t_qb = TFPQ_uncorrected_t_qb;
    ws_plots.TFPQ_uncorrected_t_q1 = TFPQ_uncorrected_t_q1;
    ws_plots.M_norm_t_qb = M_norm_t_qb;
    ws_plots.M_norm_t_q1 = M_norm_t_q1;

    % sd mpk
    ws_plots.sdmpk_t_q1 = sdmpk_t_q1;
    ws_plots.sdmpk_t_qb = sdmpk_t_qb;


    % Thresholds
    ws_plots.out_qb = out_qb;
    ws_plots.out_q1 = out_q1;
    ws_plots.out_qb_trade = out_qb_trade;
    ws_plots.out_q1_trade = out_q1_trade;

    % I*(I<0)/K (aggregate I / aggregate k)
    ws_plots.IK_neg_q1 = inc_Q_neg_t_q1./k_t_q1;
    ws_plots.IK_neg_qb = inc_Q_neg_t_qb./k_t_qb;
    ws_plots.IK_pos_q1 = inc_Q_plus_t_q1./k_t_q1;
    ws_plots.IK_pos_qb = inc_Q_plus_t_qb./k_t_qb;

    % I(I==0)
    ws_plots.inc_Q_0_t_qb = inc_Q_0_t_qb;
    ws_plots.inc_Q_0_t_q1 = inc_Q_0_t_q1;
    
    % Inaction 10%
    ws_plots.inc_Q_pos_10p_qb = inc_Q_pos_10p_qb;
    ws_plots.inc_Q_neg_10p_qb = inc_Q_neg_10p_qb;
    ws_plots.inc_Q_0_10p_qb = inc_Q_0_10p_qb;

    % I*(I<0)/K (intensive margin I / aggregate k)
    ws_plots.IK_neg_intensive_qb = (inv_neg_intensive_t_qb.*Mstay_t_qb)./k_t_qb;
    ws_plots.IK_pos_intensive_qb = (inv_pos_intensive_t_qb.*Mstay_t_qb)./k_t_qb;
    ws_plots.IK_neg_intensive_q1 = (inv_neg_intensive_t_q1.*Mstay_t_q1)./k_t_q1;
    ws_plots.IK_pos_intensive_q1 = (inv_pos_intensive_t_q1.*Mstay_t_q1)./k_t_q1;

    save('ws_plots','ws_plots','gr_*','out_qb_trade','out_q1_trade','stats_out','tab_*');

catch

    disp('Only exporting q<Q')

    % prices
    ws_plots.P_plot_qb = P_plot_qb;

    % cons
    ws_plots.C_t_qb = C_t_qb;

    % capital
    ws_plots.k_t_qb = k_t_qb;

    % Labor manu
    ws_plots.l_t_qb = l_t_qb;

    % startup rate
    ws_plots.Mentry_t_qb = Mentry_t_qb;

    % exit rate
    ws_plots.Mexit_t_qb = Mexit_t_qb;

    % Mass of active firms
    ws_plots.MM_t_qb = MM_t_qb;

    % Mass of continuing firms
    ws_plots.MM_t_qb = Mstay_t_qb;

    % Investment
    ws_plots.inc_Q_t_qb = inc_Q_t_qb;

    % Import penetration
    ws_plots.imp_t_qb = P_Cf_t_qb(:)./(P_Cf_t_qb(:)+P_Cd_t_qb(:));

    % Imports relative to domestic output
    ws_plots.rel_Cf_Cd_qb = P_Cf_t_qb(:)./P_Cd_t_qb(:);

    % TFP (raw average)
    ws_plots.avg_raw_s_t_qb = avg_raw_s_t_qb;

    % TFP (PY weighted average)
    ws_plots.avg_wgt_s_t_qb = avg_wgt_s_t_qb;

    % TFP (Y weighted average)
    ws_plots.avg_wgt_y_s_t_qb = avg_wgt_y_s_t_qb;

    % TFP (K weighted average)
    ws_plots.avg_wgt_k_s_t_qb = avg_wgt_k_s_t_qb;

    % Measured TFPQ
    ws_plots.tfpq_t_qb = tfpq_t_qb;

    % Measured TFPR
    ws_plots.tfpr_t_qb = tfp_t_qb;
    
    % True physical TFPQ
    ws_plots.TFPQ_t_qb = TFPQ_t_qb;

    % sd mpk
    ws_plots.sdmpk_t_qb = sdmpk_t_qb;

    % Thresholds
    ws_plots.out_qb = out_qb;
    ws_plots.out_qb_trade = out_qb_trade;

    % I*(I<0)/K (aggregate I / aggregate k)
    ws_plots.IK_neg_qb = inc_Q_neg_t_qb./k_t_qb;
    ws_plots.IK_pos_qb = inc_Q_plus_t_qb./k_t_qb;

    % I(I==0)
    ws_plots.inc_Q_0_t_qb = inc_Q_0_t_qb;

    % I*(I<0)/K (intensive margin I / aggregate k)
    ws_plots.IK_neg_intensive_qb = (inv_neg_intensive_t_qb.*Mstay_t_qb)./k_t_qb;
    ws_plots.IK_pos_intensive_qb = (inv_pos_intensive_t_qb.*Mstay_t_qb)./k_t_qb;

    save('ws_plots','ws_plots','gr_*','out_qb_trade');

end

% return to current path
cd(pathnow)






