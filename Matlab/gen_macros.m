% load i/k vs mrpk regressions
    % baseline
load('mat_files_compiled/baseline/mat_files/results/event_plots','beta_hat','Nsample','STATS')
beta_hat_baseline = beta_hat;
Nsample_baseline = Nsample;
se_baseline = STATS.SE;
clear beta_hat Nsample STATS
    % frictionless
load('mat_files_compiled/baseline/mat_files/results/event_plots_q1','beta_hat','Nsample','STATS')
beta_hat_frictionless = beta_hat;
Nsample_frictionless = Nsample;
se_frictionless = STATS.SE;
clear beta_hat Nsample STATS

% load event study
load('mat_files_compiled/baseline/mat_files/results/event_plots')

% load decomposition analysis
    % baseline
load('mat_files_compiled/Recreated_workspace/mat_files/transition_recreate_baseline', 'TFPQ_counterfactual_6_t', 'TFPQ_t','sd_mpk_6_t','sd_mpk_t')
TFPQ_counterfactual_6_t_baseline = TFPQ_counterfactual_6_t;
TFPQ_t_baseline = TFPQ_t;
sd_mpk_6_t_baseline = sd_mpk_6_t;
sd_mpk_t_baseline = sd_mpk_t;
clear TFPQ_counterfactual_6_t TFPQ_t sd_mpk_6_t sd_mpk_t
    % no eta
load('mat_files_compiled/Recreated_workspace/mat_files/transition_recreate_no_eta', 'TFPQ_counterfactual_6_t', 'TFPQ_t','sd_mpk_6_t','sd_mpk_t')
TFPQ_counterfactual_6_t_no_eta = TFPQ_counterfactual_6_t;
TFPQ_t_no_eta = TFPQ_t;
sd_mpk_6_t_no_eta = sd_mpk_6_t;
sd_mpk_t_no_eta = sd_mpk_t;
clear TFPQ_counterfactual_6_t TFPQ_t sd_mpk_6_t sd_mpk_t

% make directory to store macros file
mkdir('paper_numbers')

%% Construct table of macros %%%%%%
FID = fopen('table_macros_model.tex','w');
% page 21 (% capital resale loss)
fprintf(FID,' %% Percent lambda \n');
fprintf(FID, ['\\newcommand{\\perLambda}{' num2str(round(100*(par2_baseline.lam),0)) '\\%%' '} \n']);
fprintf(FID,' %% Percent zeta \n');
fprintf(FID, ['\\newcommand{\\perZeta}{' num2str(round(100*(par2_baseline.zet),0)) '\\%%' '} \n']);
% page 23 (mpk elasticity)
fprintf(FID,' %% i/k elasticity \n');
fprintf(FID, ['\\newcommand{\\ikVmpk}{' num2str(round(beta_hat_baseline(2),3),'%0.3f') '} \n']);
fprintf(FID,' %% i/k elasticity 2 decimal places \n');
fprintf(FID, ['\\newcommand{\\ikVmpkTwoDP}{' num2str(round(beta_hat_baseline(2),2),'%0.2f') '} \n']);
    % frictionless
fprintf(FID,' %% i/k elasticity, frictionless \n');
fprintf(FID, ['\\newcommand{\\ikVmpkF}{' num2str(round(beta_hat_frictionless(2),3)) '} \n']);
    % as factor of frictionless
fprintf(FID,' %% i/k elasticity, frictionless / baseline \n');
fprintf(FID, ['\\newcommand{\\ikVmpkFactor}{' num2words(round(beta_hat_frictionless(2)/beta_hat_baseline(2),0)) '} \n']);
    % sample size, baseline
fprintf(FID,' %% SE of i/k elasticity reg \n');
fprintf(FID, ['\\newcommand{\\SEikreg}{' num2str(round(se_baseline(2),3)) '} \n']);
    % SE, frictionless
fprintf(FID,' %% SE of i/k elasticity reg, frictionless \n');
fprintf(FID, ['\\newcommand{\\SEikregF}{' num2str(round(se_frictionless(2),3)) '} \n']);
% sample size, baseline
fprintf(FID,' %% sample size of i/k elasticity reg \n');
fprintf(FID, ['\\newcommand{\\samplesizeikreg}{' num2str(Nsample_baseline) '} \n']);
    % sample size, frictionless
fprintf(FID,' %% sample size of i/k elasticity reg, frictionless \n');
fprintf(FID, ['\\newcommand{\\samplesizeikregF}{' num2str(Nsample_frictionless) '} \n']);
% lumpiness
    % 20%
fprintf(FID,' %% lumpiness 20%% \n');
fprintf(FID, ['\\newcommand{\\lumpyCHp}{' num2str(round(1-output_moment_ss.avg.table_ik_moments(8),2),'%0.2f') '} \n']);
    % 49%
fprintf(FID,' %% lumpiness 49%% \n');
fprintf(FID, ['\\newcommand{\\lumpyBaselinep}{' num2str(round(1-output_moment_ss.avg.table_ik_moments(9),2),'%0.2f') '} \n']);
    % <-20%
fprintf(FID,' %% lumpiness 20%% \n');
fprintf(FID, ['\\newcommand{\\lumpynegCHp}{' num2str(round(1-output_moment_ss.stats.frac_ik_less_20,2),'%0.2f') '} \n']);
    % <-49%
fprintf(FID,' %% lumpiness 49%% \n');
fprintf(FID, ['\\newcommand{\\lumpynegBaselinep}{' num2str(round(1-output_moment_ss.stats.frac_ik_less_49,2),'%0.2f') '} \n']);
% page 24 (mpk dispersion)
fprintf(FID,' %% sd mpk \n');
fprintf(FID, ['\\newcommand{\\sdmpk}{' num2str(round(output_moment_ss.stats.sd_mpk,2)) '} \n']);
    % frictionless
fprintf(FID,' %% sd mpk, frictionless \n');
fprintf(FID, ['\\newcommand{\\sdmpkF}{' num2str(round(output_moment_ss_q1.stats.sd_mpk,2)) '} \n']);
% page 23 (inaction region changes)
fprintf(FID,' %% percentage point increase in inaction \n');
fprintf(FID, ['\\newcommand{\\inactpp}{' num2str(round(100*(ws_plots_baseline.inc_Q_0_t_qb(2)-ws_plots_baseline.inc_Q_0_t_qb(1)),1),1) '} \n']);
% page 23 (inaction region positive changes)
fprintf(FID,' %% percentage point increase in positive inaction \n');
fprintf(FID, ['\\newcommand{\\inactPospp}{' num2str(-1*round(100*(ws_plots_baseline.inc_Q_pos_10p_qb(2)-ws_plots_baseline.inc_Q_pos_10p_qb(1)),1),1) '} \n']);
% page 31 table 4
fprintf(FID,' %% LR change in C \n');
fprintf(FID, ['\\newcommand{\\delC}{' num2str(round(100*(output_moment_trade.avg.C./output_moment_ss.avg.C-1),2),'%.2f') '\\%%' '} \n']);
% page 33
fprintf(FID,' %% Baseline SR TFP change \n');
fprintf(FID, ['\\newcommand{\\TFPdBaseline}{' num2str(abs(round(100*(ws_plots_baseline.TFPQ_t_qb(3)/ws_plots_baseline.TFPQ_t_qb(1)-1),1))) '} \n']);
% page 33 footnote 42
fprintf(FID,' %% Frictionless SR TFP change \n');
fprintf(FID, ['\\newcommand{\\TFPdFrictionless}{' num2str(abs(round(100*(ws_plots_baseline.TFPQ_t_q1(3)/ws_plots_baseline.TFPQ_t_q1(1)-1),1))) '} \n']);
% page 33
fprintf(FID,' %% Initial TFP ratio (baseline/frictionless) \n');
fprintf(FID, ['\\newcommand{\\ssTFPratio}{' num2str(-1*round(eff_gap(1),1)) '} \n']);
% page 33
fprintf(FID,' %% SR TFP ratio (baseline/frictionless) \n');
fprintf(FID, ['\\newcommand{\\SRTFPratio}{' num2str(-1*round(eff_gap(3),1)) '} \n']);
% page 33 footnote 42
fprintf(FID,' %% SR TFP ratio change (baseline/frictionless) \n');
fprintf(FID, ['\\newcommand{\\effgap}{' num2str(-1*round(delta_eff_gap,1)) '} \n']);
% page 38
fprintf(FID,' %% LR TFP gains \n');
fprintf(FID, ['\\newcommand{\\tfpgainsLR}{' num2str(round(100*(output_moment_trade.avg.TFPQ/output_moment_ss.avg.TFPQ-1),1)) '} \n']);
% page 33
fprintf(FID,' %% Baseline uncorrected SR TFP change \n');
fprintf(FID, ['\\newcommand{\\TFPUdBaseline}{' num2str(abs(round(100*(ws_plots_baseline.TFPQ_uncorrected_t_qb(3)/ws_plots_baseline.TFPQ_uncorrected_t_qb(1)-1),1))) '} \n']);
% page 33 footnote 42
fprintf(FID,' %% Frictionless uncorrected SR TFP change \n');
fprintf(FID, ['\\newcommand{\\TFPUdFrictionless}{' num2str(round(eff_gap_uncorrected(1),1)) '} \n']);
% page 33 (thereabouts)
fprintf(FID,' %% Initial uncorrected TFP ratio (baseline/frictionless) \n');
fprintf(FID, ['\\newcommand{\\ssTFPratioU}{' num2str(-1*round(eff_gap_uncorrected(1),0)) '} \n']);
% page 33 (thereabouts)
fprintf(FID,' %% SR uncorrected TFP ratio (baseline/frictionless) \n');
fprintf(FID, ['\\newcommand{\\SRTFPratioU}{' num2str(-1*round(eff_gap_uncorrected(3),0)) '} \n']);% page 33 footnote 42
fprintf(FID,' %% SR uncorrected TFP ratio change (baseline/frictionless) \n');
fprintf(FID, ['\\newcommand{\\effgapU}{' num2str(-1*round(delta_eff_gap_uncorrected,1)) '} \n']);
% page 38
TFPQ_uncorrected_ss = output_moment_ss.avg.Cd/(output_moment_ss.avg.Kd^par2_baseline.alph * output_moment_ss.avg.Ld.^(1-par2_baseline.alph));
TFPQ_uncorrected_trade = output_moment_trade.avg.Cd/(output_moment_trade.avg.Kd^par2_baseline.alph * output_moment_trade.avg.Ld.^(1-par2_baseline.alph));
fprintf(FID,' %% LR uncorrected TFP ratio \n');
fprintf(FID, ['\\newcommand{\\tfpgainsLRU}{' num2str(-1*round(100*(TFPQ_uncorrected_trade/TFPQ_uncorrected_ss-1),1)) '} \n']);
% page 38 (TFP change decomposition, baseline)
fprintf(FID,' %% Baseline SR TFP change: Decomposition kpr \n');
fprintf(FID, ['\\newcommand{\\TFPdBaselineKpr}{' num2str(-1*round(100*(TFPQ_counterfactual_6_t_baseline(2)/TFPQ_counterfactual_6_t_baseline(1)-1),1)) '} \n']);
    % Extensive margin, In text: contributes positively
fprintf(FID, ['\\newcommand{\\TFPdBaselineExt}{' num2str(round(100*( (TFPQ_t_baseline(2)-TFPQ_counterfactual_6_t_baseline(2))/TFPQ_counterfactual_6_t_baseline(1)),1)) '} \n']);
    % dispersion of MRPK
fprintf(FID,' %% Baseline SR dispersion MRPK change: Decomposition kpr \n');
fprintf(FID, ['\\newcommand{\\TFPdBaselineMRPK}{' num2str(round(100*(sd_mpk_6_t_baseline(2)/sd_mpk_6_t_baseline(1)-1),1)) '} \n']);
% page 38 (SR TFP change, no eta model)
fprintf(FID,' %% Baseline SR TFP change, no eta \n');
fprintf(FID, ['\\newcommand{\\TFPdNoEta}{' num2str(abs(round(100*(ws_plots_no_eta.TFPQ_t_qb(3)/ws_plots_no_eta.TFPQ_t_qb(1)-1),1))) '} \n']);
fprintf(FID, ['\\newcommand{\\TFPdNoEtaFrictionless}{' num2str(abs(round(100*(ws_plots_no_eta.TFPQ_t_q1(3)/ws_plots_no_eta.TFPQ_t_q1(1)-1),1))) '} \n']);
% page 38 (TFP change decomposition, no eta)
fprintf(FID,' %% Baseline SR TFP change, no eta: Decomposition kpr \n');
fprintf(FID, ['\\newcommand{\\TFPdNoEtaKpr}{' num2str(-1*round(100*(TFPQ_counterfactual_6_t_no_eta(2)/TFPQ_counterfactual_6_t_no_eta(1)-1),1)) '} \n']);
    % Extensive margin, In text: contributes positively
fprintf(FID, ['\\newcommand{\\TFPdNoEtaExt}{' num2str(round(100*( (TFPQ_t_no_eta(2)-TFPQ_counterfactual_6_t_no_eta(2))/TFPQ_counterfactual_6_t_no_eta(1)),2)) '} \n']);
% page 33 (change in TFP ratio, no eta)
fprintf(FID,' %% SR TFP ratio change (baseline/frictionless), no eta \n');
fprintf(FID, ['\\newcommand{\\effgapNoEta}{' num2str(-1*round(delta_eff_gap_fc,1)) '} \n']);
% page 39 (welfare)
fprintf(FID,' %% welfare gains incl transition \n');
fprintf(FID, ['\\newcommand{\\cev}{' num2str(round(100*(ws_plots_baseline.tab_welfare(2,2)),2)) '\\%%' '} \n']);
% page 45
fprintf(FID,' %% SR TFP ratio, PE \n');
fprintf(FID, ['\\newcommand{\\effgapSRPE}{' num2str(-1*round(100*(eff_gap_PE(3)),1)) '\\%%' '} \n']);
% page 46
fprintf(FID,' %% SR TFP ratio, no convex \n');
fprintf(FID, ['\\newcommand{\\effgapSRNC}{' num2str(-1*round(100*(eff_gap_nc(3)),1)) '\\%%' '} \n']);
% page 46
fprintf(FID,' %% SR TFP ratio, eta1=eta2=0 \n');
fprintf(FID, ['\\newcommand{\\effgapSRFC}{' num2str(-1*round(100*(eff_gap_fc(3)),1)) '\\%%' '} \n']);
% re-used macros
fprintf(FID,' %% parameter value: high sigma \n');
fprintf(FID, ['\\newcommand{\\highsig}{' num2str(round(par2_high_sigma.sigma_us,3),3) '} \n']);
% fprintf(FID,' %% SR TFP gains \n');
% fprintf(FID, ['\\newcommand{\\tfpgainsSR}{' num2str(round(-100*((ws_plots_baseline.TFPQ_t_qb(3)/ws_plots_baseline.TFPQ_t_qb(2)-1)),2)) '\\%%' '} \n']);
fprintf(FID,' %% LR avg tfp gains \n');
fprintf(FID, ['\\newcommand{\\avgSgainsLR}{' num2str(round(100*((ws_plots_baseline.avg_raw_s_t_qb(end)/ws_plots_baseline.avg_raw_s_t_qb(1)-1)),1)) '\\%%' '} \n']);
fprintf(FID,' %% SR avg tfp gains \n');
fprintf(FID, ['\\newcommand{\\avgSgainsSR}{' num2str(round(-100*((ws_plots_baseline.avg_raw_s_t_qb(3)/ws_plots_baseline.avg_raw_s_t_qb(2)-1)),2)) '\\%%' '} \n']);
fprintf(FID,' %% SR exit rate change \n');
fprintf(FID, ['\\newcommand{\\exitSR}{' num2str(round(100*((exit_rate_t_qb(2)-exit_rate_t_qb(1))),1)) '} \n']);
fprintf(FID,' %% SR exit rate change, Frictionless \n');
fprintf(FID, ['\\newcommand{\\exitSRfrictionless}{' num2str(round(100*((exit_rate_t_q1(2)-exit_rate_t_q1(1))),1)) '} \n']);
fprintf(FID,' %% SR entry rate change \n');
fprintf(FID, ['\\newcommand{\\entrySR}{' num2str(round(-100*((entry_rate_t_qb(2)-entry_rate_t_qb(1))),1)) '} \n']);
fprintf(FID,' %% SR entry rate change, Frictionless \n');
fprintf(FID, ['\\newcommand{\\entrySRfrictionless}{' num2str(round(-100*((entry_rate_t_q1(2)-entry_rate_t_q1(1))),1)) '} \n']);
fprintf(FID,' %% Lumpiness moment \n');
fprintf(FID, ['\\newcommand{\\lumpiness}{' num2str(round(100*(1-output_moment_ss.stats.ik_inaction_050),2),'%.2f') '\\%%'  '} \n']);
% event study macros
fprintf(FID,' %% prop of firms ever having 1 spike \n');
fprintf(FID, ['\\newcommand{\\spikeproportionONE}{' num2str(round(pop_firms_spike1/pop_firms,3),'%.3f')   '} \n']);
fprintf(FID,' %% prop of firms having >=2 spikes \n');
fprintf(FID, ['\\newcommand{\\spikeproportionTWO}{' num2str(round(pop_firms_spike2/pop_firms,3),'%.3f')   '} \n']);
fprintf(FID,' %% Average time gap till adjustment, condn on >=2 spikes \n');
fprintf(FID, ['\\newcommand{\\timegaptTWO}{' num2str(round(mean(t_gap_spike,'omitnan'),3),'%.3f')   '} \n']);
fprintf(FID,' %% Average capital age \n');
fprintf(FID, ['\\newcommand{\\age}{' num2str(round(mean(t_age(spike2>=1),'omitnan'),3),'%.3f')   '} \n']);
fprintf(FID,' %% Average capital age, condn on only 1 spike \n');
fprintf(FID, ['\\newcommand{\\ageONEonly}{' num2str(round(mean(t_gap(spike2==1),'omitnan'),3),'%.3f')   '} \n']);
fclose(FID);

i_success = movefile('table_macros_model.tex','paper_numbers');

%%  Construct table of macros (data) %%%%%%
FID = fopen('table_macros_data.tex','w');
% page 21 (% capital resale loss)
fprintf(FID,' %% Percent lambda \n');

% SD MRPK
fprintf(FID, ['\\newcommand{\\sdmprk}{1.47} \n']);
% SD MRPK 4digits
fprintf(FID, ['\\newcommand{\\sdmprkqd}{1.42} \n']);
% MRPK Autocorrelation
fprintf(FID, ['\\newcommand{\\rhomrpk}{0.73} \n']);
% Percent negative adjustments
fprintf(FID, ['\\newcommand{\\neginv}{11\\%% } \n']);
% Transition Matrix - bottom tercile
fprintf(FID, ['\\newcommand{\\tmbottom}{82\\%%} \n']);
% Transition Matrix - top tercile
fprintf(FID, ['\\newcommand{\\tmtop}{77\\%%} \n']);
% Survival - logomega
fprintf(FID, ['\\newcommand{\\survomega}{0.29 \\; (0.02) } \n']);
% Survival - logk
fprintf(FID, ['\\newcommand{\\survk}{0.21 \\; (0.01) } \n']);
 % i/k elasticity data
fprintf(FID, ['\\newcommand{\\ikVmpkdata}{0.48 } \n']);

i_success = movefile('table_macros_data.tex','paper_numbers');
