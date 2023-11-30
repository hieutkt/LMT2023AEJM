clear
clc
close all

%% addpaths

addpath('functions');
addpath('utility functions');

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
try
%%%%%%%% steady-state %%%%%%%%%%%

% baseline economy
load('mat_files_compiled/baseline/mat_files/ss_autarky_baseline.mat','par2','output_moment')
output_moment_ss = output_moment;
par2_baseline = par2;
load('mat_files_compiled/baseline/mat_files/ss_tradeshock_baseline.mat','output_moment')
output_moment_trade = output_moment;
load('mat_files_compiled/baseline/mat_files/ss_autarky_q1_baseline.mat','output_moment')
output_moment_ss_q1 = output_moment;

% no convex
load('mat_files_compiled/no_convex/mat_files/ss_autarky_no_convex.mat','par2')
par2_no_convex = par2;

% no eta
load('mat_files_compiled/no_eta/mat_files/ss_autarky_no_eta.mat','par2')
par2_no_eta = par2;

% high sigma
load('mat_files_compiled/high_sigma/mat_files/ss_autarky_high_sigma.mat','par2')
par2_high_sigma = par2;

% fixed cost
load('mat_files_compiled/fixed_cost/mat_files/stats_out_fixedcost.mat','moments_outcome','par2')
moments_outcome_fixedcost = moments_outcome;
par2_fixedcost = par2;
clear moments_outcome par2

%%%%%%%% transitions %%%%%%%%%%%

% baseline economy
load('mat_files_compiled/baseline/mat_files/ws_plots.mat','ws_plots','stats_out','gr_qb','gr_q1')
ws_plots_baseline = ws_plots;
stats_out_baseline = stats_out;
out_qb = ws_plots_baseline.out_qb;
out_qb_trade = ws_plots_baseline.out_qb_trade;
out_q1 = ws_plots_baseline.out_q1;
out_q1_trade = ws_plots_baseline.out_q1_trade;
gr_qb_baseline = gr_qb;
gr_q1_baseline = gr_q1;

% no convex
load('mat_files_compiled/no_convex/mat_files/ws_plots.mat','ws_plots','stats_out','gr_qb')
ws_plots_no_convex = ws_plots;
stats_out_no_convex = stats_out;
out_qb_no_convex = ws_plots_no_convex.out_qb;
out_qb_no_convex_trade = ws_plots_no_convex.out_qb_trade;
gr_qb_no_convex = gr_qb;

% no eta
load('mat_files_compiled/no_eta/mat_files/ws_plots.mat','ws_plots','stats_out','gr_qb')
stats_out_no_eta = stats_out;
ws_plots_no_eta = ws_plots;
out_qb_no_eta = ws_plots_no_eta.out_qb;
out_qb_no_eta_trade = ws_plots_no_eta.out_qb_trade;
gr_qb_no_eta = gr_qb;

% PE
% load('../0_paper_convex_ik/counterfactual_cali_PE/PE_mat.mat')
load('mat_files_compiled/baseline_PE/mat_files/ws_plots.mat','ws_plots','gr_qb','gr_q1')
ws_plots_PE = ws_plots;
out_qb_PE = ws_plots_PE.out_qb_trade;
gr_qb_PE = gr_qb;
gr_q1_PE = gr_q1;

% high sigma economy
load('mat_files_compiled/high_sigma/mat_files/ws_plots.mat','ws_plots','stats_out','gr_qb')
stats_out_high_sigma = stats_out;
ws_plots_high_sigma = ws_plots;
out_qb_high_sigma = ws_plots_high_sigma.out_qb;
out_qb_high_sigma_trade = ws_plots_high_sigma.out_qb_trade;
gr_qb_high_sigma = gr_qb;

% simple figure correction
    % replace q=1 t=2 plot
ws_plots_high_sigma.TFPQ_t_qb(2) = ws_plots_high_sigma.TFPQ_t_qb(1);
ws_plots_no_convex.TFPQ_t_qb(2) = ws_plots_no_convex.TFPQ_t_qb(1);
ws_plots_baseline.TFPQ_t_qb(2) = ws_plots_baseline.TFPQ_t_qb(1);
ws_plots_no_eta.TFPQ_t_qb(2) = ws_plots_no_eta.TFPQ_t_qb(1);
ws_plots_PE.TFPQ_t_qb(2) = ws_plots_PE.TFPQ_t_qb(1);

ws_plots_high_sigma.TFPQ_t_q1(2) = ws_plots_high_sigma.TFPQ_t_q1(1);
ws_plots_no_convex.TFPQ_t_q1(2) = ws_plots_no_convex.TFPQ_t_q1(1);
ws_plots_baseline.TFPQ_t_q1(2) = ws_plots_baseline.TFPQ_t_q1(1);
ws_plots_no_eta.TFPQ_t_q1(2) = ws_plots_no_eta.TFPQ_t_q1(1);

% set plot parameters
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

% Efficiency gap calculations: No eta
    % difference in pp
TFPQ_diff_fc = ws_plots_no_eta.TFPQ_t_qb(1:T) - ws_plots_no_eta.TFPQ_t_q1(1:T);
    % difference normalized by period 1 baseline
eff_gap_fc = 100*TFPQ_diff_fc/ws_plots_no_eta.TFPQ_t_qb(1);
    % change in eff gap
delta_eff_gap_fc = eff_gap_fc(3)-eff_gap_fc(2);

%% Main text figures

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
print('figures/fig_01_baseline_inv_exit_qb','-depsc')

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
print('figures/fig_02_baseline_inv_exit_q1','-depsc')

%%% Figure 2 %%%

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
print('figures/fig_03_baseline_P','-depsc')

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
print('figures/fig_04_baseline_K','-depsc')

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
print('figures/fig_05_Y_domestic','-depsc')

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
print('figures/fig_06_baseline_M','-depsc')

%%% Figure 3 %%%

% Figure 3a: Agg TFP
figure
hold on
plot(0:T-1,100*(ws_plots_baseline.TFPQ_t_qb(1:T)/ws_plots_baseline.TFPQ_t_qb(1)-1), 'Color', xblue,'linewidth',3)
plot(0:T-1,100*(ws_plots_baseline.TFPQ_t_q1(1:T)/ws_plots_baseline.TFPQ_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('$TFP^{Adj}$($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_07_baseline_TFP','-depsc')

% Figure 3b: sd mrpk
figure
hold on
plot(0:T-1,100*(ws_plots_baseline.sdmpk_t_qb(1:T)/ws_plots_baseline.sdmpk_t_qb(1)-1), 'Color', xblue,'linewidth',3)
plot(0:T-1,100*(ws_plots_baseline.sdmpk_t_q1(1:T)/ws_plots_baseline.sdmpk_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('s.d. of MRPK ($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','northeast')
grid on
print('figures/fig_08_baseline_sd_mrpk','-depsc')

%%% Figure 4 %%%

% Figure 4a: Inaction region Before / After trade shock
s_new = (gr_qb_baseline.s_grid(1)):.1:gr_qb_baseline.s_grid(end);
k_up_new = interp1((gr_qb_baseline.s_grid),out_qb.k_up_threshold,s_new,'pchip','extrap');
k_up_new_trade = interp1((gr_qb_baseline.s_grid),out_qb_trade.k_up_threshold,s_new,'pchip','extrap');
% plot thresholds: q<Q
figure
hold on
p1 = plot(k_up_new,s_new, 'Color', xblue, 'LineWidth', 3);
p2 = plot(out_qb.k_down_threshold,(gr_qb_baseline.s_grid), 'Color', xblue, 'LineWidth', 2);
p3 = plot(k_up_new_trade,s_new, '--', 'Color', xblue, 'LineWidth', 3);
p4 = plot(out_qb_trade.k_down_threshold,(gr_qb_baseline.s_grid), '--', 'Color', xblue, 'LineWidth', 2);
hold off
xmin = 4.5;
xmax = 7.5;
ymin = 0.5;
ymax = 11.1;
ydif = ymax-ymin;
% xmin = 0;
% xmax = gr_qb.k_grid(end)-.1;
axis([xmin xmax-.1 ymin ymax])
grid on
xlabel('$k$','FontSize',24,'interpreter','latex')
ylabel('$s$','FontSize',24,'interpreter','latex')
legend({'$t=0$','$t=0$','$t=1$','$t=1$'},'location','northwest','Fontsize',20,'interpreter','latex')
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
print('figures/fig_09_baseline_inaction','-depsc')

% Figure 4b: Exit
figure
hold on
plot(gr_qb_baseline.k_grid,(out_qb.s_fx_star), 'Color', xblue, 'LineWidth', 3);
plot(gr_q1_baseline.k_grid,(out_q1.s_fx_star), 'Color', xred, 'LineWidth', 1);
plot(gr_qb_baseline.k_grid,(out_qb_trade.s_fx_star), '--', 'Color', xblue,'LineWidth', 3);
plot(gr_q1_baseline.k_grid,(out_q1_trade.s_fx_star), '--', 'Color', xred,'LineWidth', 1);
xlabel('$k$','FontSize',24,'interpreter','latex')
ylabel('$s$','FontSize',24,'interpreter','latex')
legend({'Baseline $t=0$', 'Frictionless $t=0$', 'Baseline, $t=1$', 'Frictionless, $t=1$'},'FontSize',20,'interpreter','latex')
hold off
xmin = 0;
xmax = 0.8;
ymin = 0.35;
ymax = 1.05;
axis([0 xmax-.1 ymin ymax])
grid on
print('figures/fig_10_baseline_exit','-depsc')

%%% Figure 5 %%%

% Figure 5a: Uncorrected TFP
figure
    % graphical correction
ws_plots_baseline.TFPQ_uncorrected_t_q1(2) = ws_plots_baseline.TFPQ_uncorrected_t_q1(1);
hold on
plot(0:T-1,100*(ws_plots_baseline.TFPQ_uncorrected_t_qb(1:T)/ws_plots_baseline.TFPQ_uncorrected_t_qb(1)-1), 'Color', xblue,'linewidth',3)
plot(0:T-1,100*(ws_plots_baseline.TFPQ_uncorrected_t_q1(1:T)/ws_plots_baseline.TFPQ_uncorrected_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('TFP ($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','northeast')
grid on
print('figures/fig_11_baseline_TFP_uncorrected','-depsc')

% Figure 5b: Normalizing measure
figure
hold on
plot(0:T-1,100*(1./(ws_plots_baseline.M_norm_t_qb(1:T)/ws_plots_baseline.M_norm_t_qb(1))-1), 'Color', xblue,'linewidth',3)
plot(0:T-1,100*(1./(ws_plots_baseline.M_norm_t_q1(1:T)/ws_plots_baseline.M_norm_t_q1(1))-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('Normalizing Measure ($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_12_baseline_normalizing_m','-depsc')

%%% Figure 6 %%%

% Figure 6a: Fixed R Inaction region Before / After trade shock
s_new = (gr_qb_PE.s_grid(1)):.1:gr_qb_PE.s_grid(end);
k_up_new = interp1((gr_qb_baseline.s_grid),out_qb.k_up_threshold,s_new,'pchip','extrap');
k_up_new_trade = interp1((gr_qb_baseline.s_grid),out_qb_trade.k_up_threshold,s_new,'pchip','extrap');
k_up_new_trade_PE = interp1((gr_qb_PE.s_grid),out_qb_PE.k_up_threshold,s_new,'pchip','extrap');
% plot thresholds: q<Q
figure
hold on
p1 = plot(k_up_new,s_new, 'Color', xblue, 'LineWidth', 3);
p2 = plot(out_qb.k_down_threshold,(gr_qb_PE.s_grid), 'Color', xblue, 'LineWidth', 2);
p3 = plot(k_up_new_trade,s_new, '--', 'Color', xblue, 'LineWidth', 3);
p4 = plot(out_qb_trade.k_down_threshold,(gr_qb_PE.s_grid), '--', 'Color', xblue, 'LineWidth', 2);
p5 = plot(k_up_new_trade_PE,s_new, '-.', 'Color', xred, 'LineWidth', 3);
p6 = plot(out_qb_PE.k_down_threshold,(gr_qb_PE.s_grid), '-.', 'Color', xred, 'LineWidth', 2);
hold off
xmin = 4.5;
xmax = 7.5;
ymin = 0.5;
ymax = 11.1;
ydif = ymax-ymin;
% xmin = 0;
% xmax = gr_qb.k_grid(end)-.1;
axis([xmin xmax-.1 ymin ymax])
grid on
xlabel('$k$','FontSize',24,'interpreter','latex')
ylabel('$s$','FontSize',24,'interpreter','latex')
% legend({'$t=0$','$t=0$','$t=1$','$t=1$'},'location','best','Fontsize',20,'interpreter','latex')
legend({'$t=0$','$t=0$','Baseline, $t=1$','Baseline, $t=1$','PE, $t=1$','PE, $t=1$',},'location','east','Fontsize',20,'interpreter','latex') % new legend and scales
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
print('figures/fig_13_inaction_PE','-depsc')

% Figure 6b: Fixed R Exit Thresholds
figure
hold on
plot(gr_qb_PE.k_grid,(out_qb.s_fx_star), 'Color', xblue, 'LineWidth', 3);
plot(gr_qb_PE.k_grid,(out_qb_trade.s_fx_star), '--', 'Color', xblue,'LineWidth', 3);
plot(gr_qb_PE.k_grid,(out_qb_PE.s_fx_star), '-.', 'Color', xred, 'LineWidth', 3);
xlabel('$k$','FontSize',24,'interpreter','latex')
ylabel('$s$','FontSize',24,'interpreter','latex')
legend({'$t=0$','Baseline, $t=1$', 'PE, $t=1$'},'FontSize',20,'interpreter','latex')
hold off
xmin = 0;
xmax = 0.8;
ymin = 0.35;
ymax = 1.05;
axis([0 xmax-.1 ymin ymax])
grid on
print('figures/fig_14_exit_PE','-depsc')

%% Main text tables

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

% Table 2b: transition matrix
FID = fopen('tables/table_02_transition_matrix_baseline.tex','w');
fprintf(FID,' \\begin{tabular}{cc|*{3}{c}} \n');
fprintf(FID,' \\toprule \n');
fprintf(FID,' & & \\multicolumn{3}{c}{ at $t+1$} \\\\ \n');
fprintf(FID,' & & 1 & 2 & 3 \\\\ \n');
fprintf(FID,' \\cline{1-5}  \n');
fprintf(FID,' \\multirow{3}{*}{\\makecell{Tercile at $t$}}&  \n');
fprintf(FID,[' 1  &'  num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(1,1),2),'%.2f') '&' num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(1,2),2),'%.2f') '&' num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(1,3),2),'%.2f') '\\\\  \n']);
fprintf(FID,['  &'  '&' '&' '\\\\  \n']);
fprintf(FID,[' & 2  &'  num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(2,1),2),'%.2f') '&' num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(2,2),2),'%.2f') '&' num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(2,3),2),'%.2f') '\\\\  \n']);
fprintf(FID,['  &' '&' '&' '\\\\  \n']);
fprintf(FID,[' & 3  &'  num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(3,1),2),'%.2f') '&' num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(3,2),2),'%.2f') '&' num2str(round(output_moment_ss.stats.pr_mpk_raw_norm(3,3),2),'%.2f') '\\\\  \n']);
fprintf(FID,['  &'    '&' '&' '\\\\  \n']);
fprintf(FID,' \\bottomrule \n');
fprintf(FID,' \\end{tabular} \n');
fclose(FID);

%% Appendix Figures

%%% Figure C1 (Andrea) %%%


%%% Figure C2 (Andrea) %%%


%%% Figure D2 %%%
% Import penetration
figure
plot(1:T-1,100*(ws_plots_baseline.imp_t_qb(2:T)), 'Color', xblue,'linewidth',3)
ylabel('Import Penetration (\%)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
grid on
axis([1 T-1 100*ws_plots_baseline.imp_t_qb(2) 100*ws_plots_baseline.imp_t_qb(end)])
print('figures/fig_15_baseline_imp','-depsc')

%%% Figure D3 %%
% Manufacturing labor
figure
hold on
plot(0:T-1,100*(ws_plots_baseline.l_t_qb(1:T)/ws_plots_baseline.l_t_qb(1)-1), 'Color', xblue,'linewidth',3)
plot(0:T-1,100*(ws_plots_baseline.l_t_q1(1:T)/ws_plots_baseline.l_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
ylabel('$N$ Manufacturing (\% $\Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','northeast')
grid on
axis([0 T-1 -10 0])
hold off
print('figures/fig_16_baseline_N_manu','-depsc')

%%% Figure D4 %%
% Average productivity
figure
hold on
plot(0:T-1,100*(ws_plots_baseline.avg_raw_s_t_qb(1:T)/ws_plots_baseline.avg_raw_s_t_qb(1)-1), 'Color', xblue,'linewidth',3)
plot(0:T-1,100*(ws_plots_baseline.avg_raw_s_t_q1(1:T)/ws_plots_baseline.avg_raw_s_t_q1(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('Average $s$ ($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','Frictionless'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_17_baseline_avg_tfpq','-depsc')

%%% Figure D5 %%
% TFP^Adj PE vs baseline
    % TFPQ levels
figure
plot(0:T-1,100*(ws_plots_baseline.TFPQ_t_qb(1:T)/ws_plots_baseline.TFPQ_t_qb(1)-1), 'Color', xblue,'linewidth',3)
hold on
plot(0:T-1,100*(ws_plots_PE.TFPQ_t_qb(1:T)/ws_plots_PE.TFPQ_t_qb(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('$TFP^{Adj}$($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','PE'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_18_TFPQ_baseline_v_PE','-depsc')
   
%%% Figure D6 %%
% TFP^Adj no convex vs baseline
figure
plot(0:T-1,100*(ws_plots_baseline.TFPQ_t_qb(1:T)/ws_plots_baseline.TFPQ_t_qb(1)-1), 'Color', xblue,'linewidth',3)
hold on
plot(0:T-1,100*(ws_plots_no_convex.TFPQ_t_qb(1:T)/ws_plots_no_convex.TFPQ_t_qb(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('$TFP^{Adj}$($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','No Convex Cost'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_19_TFPQ_baseline_v_no_convex','-depsc')

%%% Figure D7 %%%

% Figure D7a: Inaction region no convex costs
% Inaction region Before / After trade shock
s_new = (gr_qb_no_convex.s_grid(1)):.1:gr_qb_no_convex.s_grid(end);
k_up_new = interp1((gr_qb_no_convex.s_grid),out_qb_no_convex.k_up_threshold,s_new,'pchip','extrap');
k_up_new_trade = interp1((gr_qb_no_convex.s_grid),out_qb_no_convex_trade.k_up_threshold,s_new,'pchip','extrap');
% plot thresholds: q<Q
figure
hold on
p1 = plot(k_up_new,s_new, 'Color', xblue, 'LineWidth', 3);
p2 = plot(out_qb_no_convex.k_down_threshold,(gr_qb_no_convex.s_grid), 'Color', xblue, 'LineWidth', 2);
p3 = plot(k_up_new_trade,s_new, '--', 'Color', xblue, 'LineWidth', 3);
p4 = plot(out_qb_no_convex_trade.k_down_threshold,(gr_qb_no_convex.s_grid), '--', 'Color', xblue, 'LineWidth', 2);
hold off
xmin = 4.5;
xmax = 7.5;
ymin = 0.5;
ymax = 11.1;
ydif = ymax-ymin;
% xmin = 0;
% xmax = gr_qb.k_grid(end)-.1;
axis([xmin xmax-.1 ymin ymax])
grid on
xlabel('$k$','FontSize',24,'interpreter','latex')
ylabel('$s$','FontSize',24,'interpreter','latex')
legend({'$t=0$','$t=0$','$t=1$','$t=1$'},'location','northwest','Fontsize',20,'interpreter','latex')
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
print('figures/fig_20_inaction_no_convex','-depsc')

% Figure D7b: Exit thresholds no convex costs
figure
hold on
plot(gr_qb_no_convex.k_grid,(out_qb_no_convex.s_fx_star), 'Color', xblue, 'LineWidth', 3);
plot(gr_qb_no_convex.k_grid,(out_qb_no_convex_trade.s_fx_star), '--', 'Color', xblue,'LineWidth', 3);
xlabel('$k$','FontSize',24,'interpreter','latex')
ylabel('$s$','FontSize',24,'interpreter','latex')
legend({'$t=0$', '$t=1$'},'FontSize',20,'interpreter','latex')
hold off
xmin = 0.1;
xmax = 0.9;
ymin = 0.40;
ymax = 0.70;
axis([xmin xmax-.1 ymin ymax])
grid on
print('figures/fig_21_exit_no_convex','-depsc')

%%% Figure D8 %%
% TFP^Adj no eta (+frictionless) vs baseline
figure
plot(0:T-1,100*(ws_plots_baseline.TFPQ_t_qb(1:T)/ws_plots_baseline.TFPQ_t_qb(1)-1), 'Color', xblue,'linewidth',3)
hold on
plot(0:T-1,100*(ws_plots_no_eta.TFPQ_t_qb(1:T)/ws_plots_no_eta.TFPQ_t_qb(1)-1), '--', 'Color', xred,'linewidth',3)
plot(0:T-1,100*(ws_plots_no_eta.TFPQ_t_q1(1:T)/ws_plots_no_eta.TFPQ_t_q1(1)-1), '-x', 'Color', 'k','linewidth',3)
hold off
ylabel('$TFP^{Adj}$($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','$\eta_{1} = \eta_{2} = 0$ (with frictions)','$\eta_{1} = \eta_{2} = 0$ (frictionless)'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_22_TFPQ_baseline_v_no_eta','-depsc')

%%% Figure D9 %%
% TFP^Adj high sigma vs baseline
figure
plot(0:T-1,100*(ws_plots_baseline.TFPQ_t_qb(1:T)/ws_plots_baseline.TFPQ_t_qb(1)-1), 'Color', xblue,'linewidth',3)
hold on
plot(0:T-1,100*(ws_plots_high_sigma.TFPQ_t_qb(1:T)/ws_plots_high_sigma.TFPQ_t_qb(1)-1), '--', 'Color', xred,'linewidth',3)
hold off
ylabel('$TFP^{Adj}$($\% \Delta$ from s.s.)','FontSize',24,'interpreter','latex')
xlabel('Years','FontSize',24,'interpreter','latex')
legend({'Baseline','High $\sigma$'},'FontSize',20,'interpreter','latex','location','southeast')
grid on
print('figures/fig_23_TFPQ_baseline_v_sigma','-depsc')

%%% Figure C1 %%%
main_figC1(xblue,xred)

%%% Figure C2 %%%
main_figC2(xblue,xred)

clc

%% Appendix Tables

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

% Table D4: Parameter values with fixed adjustment cost
FID = fopen('tables/table_05_parameters_fixedcost.tex','w');
fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
fprintf(FID,'\\begin{tabular}{cccc}  \n');
fprintf(FID,'\\hline\\hline  \n');
fprintf(FID,' Parameter & Baseline & Fixed Cost & Target / Source  \\\\  \n');
fprintf(FID,' 			\\hline  \n');
fprintf(FID,['$\\rho$ &' num2str(round(par2_baseline.rho_s,3),'%0.3f') '&' num2str(round(par2_fixedcost.rho_s,3),'%0.3f') '& Autocorrelation of $\\omega$ \\\\  \n']);
fprintf(FID,[' $\\sigma$ &' num2str(round(par2_baseline.sigma_us,3),'%0.3f') '&'  num2str(round(par2_fixedcost.sigma_us,3),'%0.3f') '& Standard deviation of $\\omega$ \\\\  \n']);
fprintf(FID,[' $q/Q$ &' num2str(round(1-par2_baseline.lam,3),'%0.3f') '&' num2str(round(1-par2_fixedcost.lam,3),'%0.3f') '& frequency of negative investment \\\\  \n']);
fprintf(FID,[' $\\zeta$ &' num2str(round(par2_baseline.zet,3),'%0.3f') '&' num2str(round(par2_fixedcost.zet,3),'%0.3f') '& Slope of exit thresholds \\\\  \n']);
fprintf(FID,[' $\\gamma_0$ &' num2str(round(par2_baseline.c1,3),'%0.3f') '&' num2str(round(par2_fixedcost.c1,3),'%0.3f') '& Standard deviation of $i/k$ \\\\  \n']);
fprintf(FID,[' $\\eta_0$ &' num2str(round(par2_baseline.mu_fstay_b,3),'%0.3f') '&' num2str(round(par2_fixedcost.mu_fstay_b,3),'%0.3f') '& Exit rate \\\\  \n']);
fprintf(FID,[' $\\eta_1$ &' num2str(round(par2_baseline.eta(1),4),'%0.4f') '&' num2str(round(par2_fixedcost.eta(1),4),'%0.4f') '& Relative size at exit \\\\  \n']);
fprintf(FID,[' $\\eta_2$ &' num2str(round(sum(par2_baseline.eta),4),'%0.4f') '&' num2str(round(sum(par2_fixedcost.eta),4),'%0.4f') '& Relative productivity at exit \\\\  \n']);
fprintf(FID,[' $f_k$ &' '-- &' num2str(round((1-par2_baseline.delta)*par2_fixedcost.C_f,4),'%0.4f') '& $Pr \\left( \\vert \\frac{i}{k} \\vert > 0.49 \\right)$ \\\\  \n']);
fprintf(FID,' \\hline\\hline  \n');
fprintf(FID,' \\end{tabular}  \n');
fprintf(FID,' }  \n');
fclose(FID);

% Table D6: Frictionless model transition matrix
FID = fopen('tables/table_transition_matrix_frictions.tex','w');
fprintf(FID,' \\begin{tabular}{cc|*{3}{c}} \n');
fprintf(FID,' \\toprule \n');
fprintf(FID,' & & \\multicolumn{3}{c}{ at $t+1$} \\\\ \n');
fprintf(FID,' & & 1 & 2 & 3 \\\\ \n');
fprintf(FID,' \\cline{1-5}  \n');
fprintf(FID,' \\multirow{3}{*}{\\makecell{Tercile at $t$}}&  \n');
fprintf(FID,[' 1  &'  num2str(0.33,'%.2f') '&' num2str(0.33,'%.2f') '&' num2str(0.33,'%.2f') '\\\\  \n']);
fprintf(FID,['  &'  '&' '&' '\\\\  \n']);
fprintf(FID,[' & 2  &'  num2str(0.33,'%.2f') '&' num2str(0.33,'%.2f') '&' num2str(0.33,'%.2f') '\\\\  \n']);
fprintf(FID,['  &' '&' '&' '\\\\  \n']);
fprintf(FID,[' & 3  &'  num2str(0.33,'%.2f') '&' num2str(0.33,'%.2f') '&' num2str(0.33,'%.2f') '\\\\  \n']);
fprintf(FID,['  &'    '&' '&' '\\\\  \n']);
fprintf(FID,' \\bottomrule \n');
fprintf(FID,' \\end{tabular} \n');
fclose(FID);

% Table D7: steady-state comparisons: trade shock LR
    % compute TFP
TFPQ_uncorrected_ss = output_moment_ss.avg.Cd/(output_moment_ss.avg.Kd^par2_baseline.alph * output_moment_ss.avg.Ld.^(1-par2_baseline.alph));
TFPQ_uncorrected_trade = output_moment_trade.avg.Cd/(output_moment_trade.avg.Kd^par2_baseline.alph * output_moment_trade.avg.Ld.^(1-par2_baseline.alph));
M_norm_ss = output_moment_ss.avg.Mactive.^(1/(par2_baseline.epsi-1));
M_norm_trade = output_moment_trade.avg.Mactive.^(1/(par2_baseline.epsi-1));
TFPQ_corrected_ss = TFPQ_uncorrected_ss/M_norm_ss;
TFPQ_corrected_trade = TFPQ_uncorrected_trade/M_norm_trade;
FID = fopen('tables/table_06_ss_agg_baseline_vs_trade.tex','w');
fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
fprintf(FID,'\\centering{} \n');
fprintf(FID,'\\begin{tabular}{cc} \n');
fprintf(FID,'\\hline\\hline \n');
fprintf(FID,'Variable & $\\Delta \\%%$ \\\\ \n');
fprintf(FID,'\\hline \n');
fprintf(FID,['$C$ &' num2str(round(100*(output_moment_trade.avg.C./output_moment_ss.avg.C-1),2),'%.2f') '\\\\ \n']);
fprintf(FID,['$K$ &' num2str(round(100*(output_moment_trade.avg.Kd./output_moment_ss.avg.Kd-1),2),'%.2f') '\\\\ \n']);
fprintf(FID,['$N$ &' num2str(round(100*(output_moment_trade.avg.Ld./output_moment_ss.avg.Ld-1),2),'%.2f') '\\\\ \n']);
fprintf(FID,['$M$ &' num2str(round(100*(output_moment_trade.avg.Mactive./output_moment_ss.avg.Mactive-1),2),'%.2f') '\\\\ \n']);
% fprintf(FID,['$TFPQ$ (average) &' num2str(round(100*(output_moment_trade.avg.avg_raw_s./output_moment_ss.avg.avg_raw_s-1),2),'%.2f') '\\\\ \n']);
fprintf(FID,['$TFP^{Adj}$ &' num2str(round(100*(TFPQ_corrected_trade./TFPQ_corrected_ss-1),2),'%.2f') '\\\\ \n']);
fprintf(FID,['$TFP$ &' num2str(round(100*(TFPQ_uncorrected_trade./TFPQ_uncorrected_ss-1),2),'%.2f') '\\\\ \n']);
fprintf(FID,'\\hline\\hline \n');
fprintf(FID,'\\end{tabular} \n');
fprintf(FID,'} \n');
fclose(FID);

% Table D8: Welfare
FID = fopen('tables/table_07_ss_welfare.tex','w');
fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
fprintf(FID,'\\centering{} \n');
fprintf(FID,'\\begin{tabular}{ccc} \n');
fprintf(FID,'\\hline\\hline \n');
fprintf(FID,'Model & Steady-state  & Transition \\\\ \n');
fprintf(FID,'\\hline \n');
fprintf(FID,['Baseline &' num2str(round(100*ws_plots_baseline.tab_welfare(2,1),2),'%.2f')  '\\%%' '&' num2str(round(100*ws_plots_baseline.tab_welfare(2,2),2),'%.2f') '\\%%' '\\\\ \n']);
fprintf(FID,['Frictionless &' num2str(round(100*ws_plots_baseline.tab_welfare(1,1),2),'%.2f') '\\%%' '&' num2str(round(100*ws_plots_baseline.tab_welfare(1,2),2),'%.2f') '\\%%' '\\\\ \n']);
fprintf(FID,'\\hline\\hline \n');
fprintf(FID,'\\end{tabular} \n');
fprintf(FID,'} \n');
fclose(FID);

% Table D9: Parameter values no convex cost
FID = fopen('tables/table_08_parameters_no_convex.tex','w');
fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
fprintf(FID,'\\begin{tabular}{ccc}  \n');
fprintf(FID,'\\hline\\hline  \n');
fprintf(FID,' Parameter & Value & Target / Source  \\\\  \n');
fprintf(FID,' 			\\hline  \n');
fprintf(FID,['$\\rho$ &' num2str(round(par2_no_convex.rho_s,3),3) '& Autocorrelation of $\\omega$ \\\\  \n']);
fprintf(FID,[' $\\sigma$ &' num2str(round(par2_no_convex.sigma_us,3),3) '& Standard deviation of $\\omega$ \\\\  \n']);
fprintf(FID,[' $q/Q$ &' num2str(round(1-par2_no_convex.lam,3),3) '& frequency of negative investment \\\\  \n']);
fprintf(FID,[' $\\zeta$ &' num2str(round(par2_no_convex.zet,3),3) '& Slope of exit thresholds \\\\  \n']);
fprintf(FID,[' $\\eta_0$ &' num2str(round(par2_no_convex.mu_fstay_b,3),3) '& Exit rate \\\\  \n']);
fprintf(FID,[' $\\eta_1$ &' num2str(round(par2_no_convex.eta(1),4),5) '& Relative size at exit \\\\  \n']);
fprintf(FID,[' $\\eta_2$ &' num2str(round(sum(par2_no_convex.eta),4),5) '& Relative productivity at exit \\\\  \n']);
fprintf(FID,' \\hline\\hline  \n');
fprintf(FID,' \\end{tabular}  \n');
fprintf(FID,' }  \n');
fclose(FID);

% Table D10: Parameter values no eta
FID = fopen('tables/table_09_parameters_no_eta.tex','w');
fprintf(FID,'{\\fontshape\\scdefault\\selectfont \n');
fprintf(FID,'\\begin{tabular}{ccc}  \n');
fprintf(FID,'\\hline\\hline  \n');
fprintf(FID,' Parameter & Value & Target / Source  \\\\  \n');
fprintf(FID,' 			\\hline  \n');
fprintf(FID,['$\\rho$ &' num2str(round(par2_no_eta.rho_s,3),3) '& Autocorrelation of $\\omega$ \\\\  \n']);
fprintf(FID,[' $\\sigma$ &' num2str(round(par2_no_eta.sigma_us,3),3) '& Standard deviation of $\\omega$ \\\\  \n']);
fprintf(FID,[' $q/Q$ &' num2str(round(1-par2_no_eta.lam,3),3) '& frequency of negative investment \\\\  \n']);
fprintf(FID,[' $\\zeta$ &' num2str(round(par2_no_eta.zet,3),3) '& Slope of exit thresholds \\\\  \n']);
fprintf(FID,[' $\\gamma_0$ &' num2str(round(par2_no_eta.c1,3),3) '& Standard deviation of $i/k$ \\\\  \n']);
fprintf(FID,[' $\\eta_0$ &' num2str(round(par2_no_eta.mu_fstay_b,3),3) '& Exit rate \\\\  \n']);
fprintf(FID,' \\hline\\hline  \n');
fprintf(FID,' \\end{tabular}  \n');
fprintf(FID,' }  \n');
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
