clear
close all

%% addpaths

addpath('functions');
addpath('utility functions');

disp("[CHY_gen_figures_tables.m] Preparing the plotting evironment...")

% removed figures folder to avoid file conflict
try
    rmdir('figures','s')
    disp('figures exist, removing folder to avoid path conflict')
catch
    disp('figures does not already exist')
end
mkdir('figures')
disp(' ')

% removed tables folder to avoid file conflict
try
    rmdir('tables','s')
    disp('tables exist, removing folder to avoid path conflict')
catch
    disp('tables does not already exist')
end
mkdir('tables')
disp(' ')

%%

%%% load files %%%
disp("[CHY_gen_figures_tables.m] Loading files...")

try
    %%%%%%%% steady-state %%%%%%%%%%%

    % baseline economy
    disp("    Loading data from the autarky model...")
    load('mat_files_compiled/baseline/mat_files/ss_autarky_baseline.mat','par2','output_moment')
    output_moment_ss = output_moment;
    par2_baseline = par2;

    disp("    Loading data from the tradeshock baseline model...")
    load('mat_files_compiled/baseline/mat_files/ss_tradeshock_baseline.mat','output_moment')
    output_moment_trade = output_moment;

    disp("    Loading data from the autarky even study model...")
    load('mat_files_compiled/baseline/mat_files/ss_autarky_q1_baseline.mat','output_moment')
    output_moment_ss_q1 = output_moment;


    %%%%%%%% transitions %%%%%%%%%%%

    % baseline economy
    disp("    Loading data from ws_plots...")
    load('mat_files_compiled/baseline/mat_files/ws_plots.mat','ws_plots','stats_out','gr_qb','gr_q1')
    ws_plots_baseline = ws_plots;
    stats_out_baseline = stats_out;
    out_qb = ws_plots_baseline.out_qb;
    out_qb_trade = ws_plots_baseline.out_qb_trade;
    out_q1 = ws_plots_baseline.out_q1;
    out_q1_trade = ws_plots_baseline.out_q1_trade;
    gr_qb_baseline = gr_qb;
    gr_q1_baseline = gr_q1;


    % set plot parameters
    disp("    Setting plot parameters...")
    xblue = [0, 0.4470, 0.7410];
    xred = [0.8500, 0.3250, 0.0980];
    T = 15;

    % TFP ratio for baseline model
    tfpq_ratio = ws_plots_baseline.TFPQ_t_qb(1:T)./ws_plots_baseline.TFPQ_t_q1(1:T);

    % Efficiency gap calculations: Adjusted
    % difference in pp
    TFPQ_diff = ws_plots_baseline.TFPQ_t_qb(1:T) - ws_plots_baseline.TFPQ_t_q1(1:T);
    % difference normalized by period 1 baseline
    eff_gap = 100*TFPQ_diff/ws_plots_baseline.TFPQ_t_qb(1);
    % change in eff gap
    delta_eff_gap = eff_gap(3)-eff_gap(2);

    % Efficiency gap calculations: Unadjusted
    % difference in pp
    TFPQ_diff_uncorrected = ws_plots_baseline.TFPQ_uncorrected_t_qb(1:T) - ws_plots_baseline.TFPQ_uncorrected_t_q1(1:T);
    % difference normalized by period 1 baseline
    eff_gap_uncorrected = 100*TFPQ_diff_uncorrected/ws_plots_baseline.TFPQ_uncorrected_t_qb(1);
    % change in eff gap
    delta_eff_gap_uncorrected = eff_gap_uncorrected(3)-eff_gap_uncorrected(2);

    % extensive margins
    % exit rate for baseline model
    exit_rate_t_qb = ws_plots_baseline.Mexit_t_qb./ws_plots_baseline.MM_t_qb;
    % entry rate for baseline model
    entry_rate_t_qb = ws_plots_baseline.Mentry_t_qb./ws_plots_baseline.MM_t_qb;
    % exit rate for frictionless model
    exit_rate_t_q1 = ws_plots_baseline.Mexit_t_q1./ws_plots_baseline.MM_t_q1;
    % entry rate for frictionless model
    entry_rate_t_q1 = ws_plots_baseline.Mentry_t_q1./ws_plots_baseline.MM_t_q1;


    %% Main text figures

    disp("[CHY_gen_figures_tables.m] Replicating Figure 1 in main text...")

    %%% Figure 1 %%%

    % Figure 1a: Thresholds for inaction and exit in s.s.: q<Q
    % subset out unique k values for downsizing (because of boundary)
    i_k_down_unqiue = out_qb.k_down_threshold<max(out_qb.k_down_threshold);
    s_down_grid = gr_qb_baseline.s_grid(i_k_down_unqiue);
    k_down_grid = out_qb.k_down_threshold(i_k_down_unqiue);
    % find crossing
    Exit_s = griddedInterpolant(gr_qb_baseline.k_grid,(out_qb.s_fx_star));
    % inv_s = griddedInterpolant(out_qb.k_up_threshold,gr_qb_baseline.s_grid);
    inv_s = griddedInterpolant(out_qb.k_up_threshold,gr_qb_baseline.s_grid);
    dinv_s = griddedInterpolant(k_down_grid,s_down_grid);
    qq_inv = @(kk) -1*(Exit_s(kk)-inv_s(kk))^2;
    qq_dinv = @(kk) -1*(Exit_s(kk)-dinv_s(kk))^2;
    Exit_inv_kstar = goldenx(qq_inv,0,6);
    Exit_dinv_kstar = goldenx(qq_dinv,0,6);
    % evaluate at crossing point and interpolate to get new plots
    Exit_inv_sstar = inv_s(Exit_inv_kstar);
    Exit_dinv_sstar = dinv_s(Exit_dinv_kstar);
    s_up_new = Exit_inv_sstar:.1:gr_qb_baseline.s_grid(end)+0.1;
    % s_down_new = Exit_dinv_sstar:.1:gr_qb_baseline.s_grid(end);
    s_down_max = s_down_grid(end);
    s_down_new = Exit_dinv_sstar:.1:s_down_max;
    k_up_new = interp1((gr_qb_baseline.s_grid),out_qb.k_up_threshold,s_up_new,'pchip','extrap');
    k_down_new = interp1(s_down_grid,k_down_grid,s_down_new,'pchip','extrap');
    % plot graph
    figure
    hold on
    plot(gr_qb_baseline.k_grid,out_qb.s_fx_star, 'LineWidth', 3);
    plot(k_up_new,s_up_new,'--', 'LineWidth', 3);
    plot(k_down_new,s_down_new,'-.', 'LineWidth', 3);
    hold off
    xlabel('$k$','FontSize',24,'interpreter','latex')
    ylabel('$s$','FontSize',24,'interpreter','latex')
    ix = find(out_qb.s_fx_star>gr_qb_baseline.s_grid(2),1,'last');
    xmin = 0;
    xmax = 1.9;
    ymin = 0;
    ymax = 2;
    axis([0 xmax-.1 ymin ymax])
    grid on
    legend({'Exit', 'i>0', 'i<0'},'location','northeast','fontsize',20)
    print('CHY_figures/fig_01_baseline_inv_exit_qb','-depsc')

    % Figure 1b: Thresholds for inaction and exit in s.s.: q=Q
    Exit_s = griddedInterpolant(gr_q1.k_grid,(out_q1.s_fx_star));
    inv_s = griddedInterpolant(gr_q1.k_grid,out_q1.s_k0_threshold);
    qq_inv = @(kk) -1*(Exit_s(kk)-inv_s(kk))^2;
    Exit_inv_kstar = goldenx(qq_inv,0,6);
    % evaluate at crossing point and interpolate to get new plots
    Exit_inv_sstar = inv_s(Exit_inv_kstar);
    s_up_new = Exit_inv_sstar:.1:gr_q1.s_grid(end)+0.1;
    k_up_new = interp1(out_q1.s_k0_threshold,gr_q1.k_grid,s_up_new,'pchip','extrap');
    % plot graph
    figure
    hold on
    plot(gr_q1.k_grid,(out_q1.s_fx_star), 'LineWidth', 3);
    plot(k_up_new,s_up_new,'--', 'LineWidth', 3);
    % plot(out_q1.s_k0_threshold,gr_q1.k_grid, 'LineWidth', 3);
    hold off
    xlabel('$k$','FontSize',24,'interpreter','latex')
    ylabel('$s$','FontSize',24,'interpreter','latex')
    xmin = 0;
    xmax = 1.9;
    ymin = 0;
    ymax = 2;
    axis([0 xmax-.1 ymin ymax])
    legend({'Exit', 'i>0'},'location','northwest','fontsize',20)
    grid on
    print('CHY_figures/fig_02_baseline_inv_exit_q1','-depsc')

    %%% Figure 2 %%%
    disp("[CHY_gen_figures_tables.m] Replicating Figure 2 in main text...")

    % Figure 2a: Prices
    figure
    hold on
    plot(0:T-1,100*(ws_plots_baseline.P_plot_qb(1:T)/ws_plots_baseline.P_plot_qb(1)-1), 'Color', xblue,'linewidth',3)
    plot(0:T-1,100*(ws_plots_baseline.P_plot_q1(1:T)/ws_plots_baseline.P_plot_q1(1)-1), '--', 'Color', xred,'linewidth',3)
    ylabel('$P$ (\% $\Delta$ from s.s.)','FontSize',24,'interpreter','latex')
    xlabel('Years','FontSize',24,'interpreter','latex')
    legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','southeast')
    grid on
    axis([0 T-1 -7 0])
    hold off
    print('CHY_figures/fig_03_baseline_P','-depsc')

    % Figure 2b: Capital
    figure
    hold on
    plot(0:T-1,100*(ws_plots_baseline.k_t_qb(1:T)/ws_plots_baseline.k_t_qb(1)-1), 'Color', xblue,'linewidth',3)
    plot(0:T-1,100*(ws_plots_baseline.k_t_q1(1:T)/ws_plots_baseline.k_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
    ylabel('$K$ (\% $\Delta$ from s.s.)','FontSize',24,'interpreter','latex')
    xlabel('Years','FontSize',24,'interpreter','latex')
    legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','northeast')
    grid on
    axis([0 T-1 -12 0])
    hold off
    print('CHY_figures/fig_04_baseline_K','-depsc')

    % Figure 2c: Domestic production
    figure
    plot(0:T-1,100*(ws_plots_baseline.Cd_t_qb(1:T)/ws_plots_baseline.Cd_t_qb(1)-1), 'Color', xblue,'linewidth',3)
    hold on
    plot(0:T-1,100*(ws_plots_baseline.Cd_t_q1(1:T)/ws_plots_baseline.Cd_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
    hold off
    ylabel('$Y^{D}$($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
    xlabel('Years','FontSize',24,'interpreter','latex')
    legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','northeast')
    grid on
    print('CHY_figures/fig_05_Y_domestic','-depsc')

    % Figure 2d: Measure of active firms
    figure
    hold on
    plot(0:T-1,100*(ws_plots_baseline.MM_t_qb(1:T)/ws_plots_baseline.MM_t_qb(1)-1), 'Color', xblue,'linewidth',3)
    plot(0:T-1,100*(ws_plots_baseline.MM_t_q1(1:T)/ws_plots_baseline.MM_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
    ylabel('$M^{D}$ (\% $\Delta$ from s.s.)','FontSize',24,'interpreter','latex')
    xlabel('Years','FontSize',24,'interpreter','latex')
    legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','northeast')
    grid on
    axis([0 T-1 -10 0])
    hold off
    print('CHY_figures/fig_06_baseline_M','-depsc')


    %% Main text tables
    disp("[CHY_gen_figures_tables.m] Replicating Table 1 in main text...")

    % Table 1: Parameter values
    FID = fopen('tables/table_01_parameters_baseline.tex','w');
    fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
    fprintf(FID,'\\begin{tabular}{ccc}  \n');
    fprintf(FID,'\\hline\\hline  \n');
    fprintf(FID,' Parameter & Value & Target / Source  \\\\  \n');
    fprintf(FID,' 			\\hline  \n');
    fprintf(FID,['$\\beta$ &' num2str(round(par2_baseline.betaa,3),3) '& Literature \\\\  \n']);
    fprintf(FID,['$\\chi$ &' num2str(round(par2_baseline.chi,3),3) '& Literature \\\\  \n']);
    fprintf(FID,['$\\epsilon$ &' num2str(round(par2_baseline.epsi,3),3) '& Literature \\\\  \n']);
    fprintf(FID,' 			\\hline  \n');
    fprintf(FID,['$\\alpha$ &' num2str(round(par2_baseline.alph,3),3) '& Capital Share \\\\  \n']);
    fprintf(FID,['$\\delta$ &' num2str(round(par2_baseline.delta,3),3) '& Depreciation Rate \\\\  \n']);
    fprintf(FID,' 			\\hline  \n');
    fprintf(FID,['$\\rho$ &' num2str(round(par2_baseline.rho_s,3),3) '& Autocorrelation of $\\omega$ \\\\  \n']);
    fprintf(FID,[' $\\sigma$ &' num2str(round(par2_baseline.sigma_us,3),3) '& Standard deviation of $\\omega$ \\\\  \n']);
    fprintf(FID,[' $q/Q$ &' num2str(round(1-par2_baseline.lam,3),3) '& frequency of negative investment \\\\  \n']);
    fprintf(FID,[' $\\zeta$ &' num2str(round(par2_baseline.zet,3),3) '& Slope of exit thresholds \\\\  \n']);
    fprintf(FID,[' $\\gamma_0$ &' num2str(round(par2_baseline.c1,3),3) '& Standard deviation of $i/k$ \\\\  \n']);
    fprintf(FID,[' $\\eta_0$ &' num2str(round(par2_baseline.mu_fstay_b,3),3) '& Exit rate \\\\  \n']);
    fprintf(FID,[' $\\eta_1$ &' num2str(round(par2_baseline.eta(1),4),5) '& Relative size at exit \\\\  \n']);
    fprintf(FID,[' $\\eta_2$ &' num2str(round(sum(par2_baseline.eta),4),5) '& Relative productivity at exit \\\\  \n']);
    fprintf(FID,' \\hline\\hline  \n');
    fprintf(FID,' \\end{tabular}  \n');
    fprintf(FID,' }  \n');
    fclose(FID);


    %% Appendix Tables
    disp("[CHY_gen_figures_tables.m] Replicating Table D1 in main text...")

    % Table D1: Calibration targets baseline
    FID = fopen('tables/table_03_calibration_baseline.tex','w');
    fprintf(FID,'\\begin{table}[H] \n');
    fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
    fprintf(FID,'		\\begin{center} \n');
    fprintf(FID,'			\\begin{tabular}{ccc} \n');
    fprintf(FID,'				\\hline\\hline \n');
    fprintf(FID,'				Moments	& Data & Model \\\\ \n');
    fprintf(FID,'				\\hline \n');
    fprintf(FID,['				Freq of negative investment & 0.108 &' num2str(round(stats_out_baseline.moments_outcome(1),3),3) '\\\\ \n']);
    fprintf(FID,['				Slope of exit thresholds & 0.754 &' num2str(round(stats_out_baseline.moments_outcome(2),3),3) '\\\\ \n']);
    fprintf(FID,['				Autocorrelation of $\\omega$ & 0.742 &' num2str(round(stats_out_baseline.moments_outcome(3),3),3) '\\\\ \n']);
    fprintf(FID,['				Unconditional std dev of $\\omega$ & 0.848 &' num2str(round(stats_out_baseline.moments_outcome(4),3),3) '\\\\ \n']);
    fprintf(FID,['				Exit rate & 0.184 &' num2str(round(stats_out_baseline.moments_outcome(5),3),3) '\\\\ \n']);
    fprintf(FID,['				Relative size at exit & 0.345 &' num2str(round(stats_out_baseline.moments_outcome(6),3),3) '\\\\ \n']);
    fprintf(FID,['				Relative productivity at exit & 0.757 &' num2str(round(stats_out_baseline.moments_outcome(7),3),3)  '\\\\ \n']);
    fprintf(FID,['				Dispersion of $\\frac{i}{k}$ & 0.828 &' num2str(round(stats_out_baseline.moments_outcome(8),3),3)  '\\\\ \n']);
    fprintf(FID,'				\\hline\\hline \n');
    fprintf(FID,'			\\end{tabular} \n');
    fprintf(FID,'	\\end{center}} \n');
    % fprintf(FID,'	\\caption{Model Fit.} \n');
    fprintf(FID,'	\\label{tab:moments-baseline} \n');
    fprintf(FID,'\\end{table} \n');
    fclose(FID);


    % Table D2: Steady-state aggregates
    disp("[CHY_gen_figures_tables.m] Replicating Table D2 in main text...")
    % compute s-s TFP
    TFPQ_uncorrected_ss = output_moment_ss.avg.Cd/(output_moment_ss.avg.Kd^par2_baseline.alph * output_moment_ss.avg.Ld.^(1-par2_baseline.alph));
    TFPQ_uncorrected_ss_q1 = output_moment_ss_q1.avg.C/(output_moment_ss_q1.avg.Kd^par2_baseline.alph * output_moment_ss_q1.avg.Ld.^(1-par2_baseline.alph));
    M_norm_ss = output_moment_ss.avg.Mactive.^(1/(par2_baseline.epsi-1));
    TFPQ_corrected_ss = TFPQ_uncorrected_ss/M_norm_ss;
    M_norm_ss_q1 = output_moment_ss_q1.avg.Mactive.^(1/(par2_baseline.epsi-1));
    TFPQ_corrected_ss_q1 = TFPQ_uncorrected_ss_q1/M_norm_ss_q1;
    FID = fopen('tables/table_04_ss_agg_baseline_vs_frictionless.tex','w');
    fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
    fprintf(FID,'\\centering{} \n');
    fprintf(FID,'\\begin{tabular}{ccc} \n');
    fprintf(FID,'\\hline\\hline \n');
    fprintf(FID,'Variable & BASELINE  & FRICTIONLESS \\\\ \n');
    fprintf(FID,'\\hline \n');
    fprintf(FID,['$C$ &' num2str(round(output_moment_ss.avg.C,2),'%.2f') '&' num2str(round(output_moment_ss_q1.avg.C,2),'%.2f') '\\\\ \n']);
    fprintf(FID,['$K$ &' num2str(round(output_moment_ss.avg.Kd,2),'%.2f') '&' num2str(round(output_moment_ss_q1.avg.Kd,2),'%.2f') '\\\\ \n']);
    fprintf(FID,['$N$ &' num2str(round(output_moment_ss.avg.L,2),'%.2f') '&' num2str(round(output_moment_ss_q1.avg.L,2),'%.2f') '\\\\ \n']);
    fprintf(FID,['$M$ &' num2str(round(output_moment_ss.avg.Mactive,2),'%.2f') '&' num2str(round(output_moment_ss_q1.avg.Mactive,2),'%.2f') '\\\\ \n']);
    fprintf(FID,['$TFP$ (average) &' num2str(round(output_moment_ss.avg.avg_raw_s,2),'%.2f') '&' num2str(round(output_moment_ss_q1.avg.avg_raw_s,2),'%.2f') '\\\\ \n']);
    fprintf(FID,['$TFP$ &' num2str(round(TFPQ_uncorrected_ss,2),'%.2f') '&' num2str(round(TFPQ_uncorrected_ss_q1,2),'%.2f') '\\\\ \n']);
    fprintf(FID,['$TFP^{Adj}$ &' num2str(round(TFPQ_corrected_ss,2),'%.2f') '&' num2str(round(TFPQ_corrected_ss_q1,2),'%.2f') '\\\\ \n']);
    fprintf(FID,'\\hline\\hline \n');
    fprintf(FID,'\\end{tabular} \n');
    fprintf(FID,'} \n');
    fclose(FID);


    %% Compute efficiency gap measures (for macros)

    % PE
    tfpq_ratio_PE = ws_plots_PE.TFPQ_t_qb(1:T)./ws_plots_PE.TFPQ_t_q1(1:T);
    TFPQ_diff_PE = ws_plots_PE.TFPQ_t_qb(1:T)/ws_plots_PE.TFPQ_t_qb(1) - ws_plots_PE.TFPQ_t_q1(1:T)/ws_plots_PE.TFPQ_t_q1(1);
    eff_gap_PE = TFPQ_diff_PE/tfpq_ratio_PE(1);

    % No convex
    tfpq_ratio_nc = ws_plots_no_convex.TFPQ_t_qb(1:T)./ws_plots_no_convex.TFPQ_t_q1(1:T);
    TFPQ_diff_nc = ws_plots_no_convex.TFPQ_t_qb(1:T)/ws_plots_no_convex.TFPQ_t_qb(1) - ws_plots_no_convex.TFPQ_t_q1(1:T)/ws_plots_no_convex.TFPQ_t_q1(1);
    eff_gap_nc = TFPQ_diff_nc/tfpq_ratio_nc(1);

catch
    disp('Error! Either mat_files_compiled not found, or certain workspaces are missing.')
    disp(' ')
    disp('Please run master.m from the parent folder to generate all required workspaces')
end
