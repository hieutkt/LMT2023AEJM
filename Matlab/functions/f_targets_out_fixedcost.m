% f_targets_out_fixedcost.m
% Compute and export s-s moments for model with fixed costs
%
% Targets:
%   1. Pr(I<0)
%   2. Exit rate (hybrid)
%   3. K ratio - exits vs all
%   4. TFP ratio - exits vs all
%   5. Slope
%   6. sd ik
%   7. rho tfpr
%   8. sd tfpr
%   9. Pr(|i/k|>x), x=0.49

function [qloss, xx, par_in] = f_targets_out_fixedcost(x0,targets,ubx,normx,wgt)

% calibrate CH model
par_in.model_fixed_cost = 'ch';

% baseline fixed parameters
    % to save output or not?
par_in.save = 0;
    % calibrating or solving model?
par_in.cali = 0;
    % fixed R (PE) or vary R (GE)?
par_in.fixedR = 0;
par_in.Me = 5;
par_in.Ns = 31;
par_in.Nk = 10000;
par_in.new_kmax = 8.0526;
par_in.new_Pnow = 1/3; % fixed P (Clementi+Palazzo calibration)
par_in.npow = 1;
par_in.kapp = 1;
par_in.s_sh = 1;
par_in.mu_fstay_a = 0;
par_in.mu_fenter_a = 0;

x0 = x0(:).*ubx(:)/normx; % rescale
par_in.lam          = x0(1);
par_in.zet          = x0(2);
par_in.c1           = x0(3);
par_in.eta(1)       = x0(4);
par_in.eta(2)       = x0(5);
par_in.mu_fstay_b   = (par_in.mu_fstay_a+x0(6));
par_in.mu_fenter_b  = (par_in.mu_fenter_a+x0(6));
par_in.rho_s        = x0(7);
par_in.sigma_us     = x0(8);
par_in.C_f          = x0(9);
    % impose identical convex cost for intensive / extensive margins
par_in.c2 = par_in.c1;

% solve model
moments_out = f_PE_ss_autarky_convex(par_in);

%%% Export targeted moments
xx = zeros(size(targets));
xx(1) = moments_out.moments.I_Gross_frac_intensive_neg;
xx(2) = moments_out.avg.Mexit/moments_out.avg.Mactive;
xx(3) = moments_out.moments.avg_k_exit/moments_out.moments.avg_k_all;
xx(4) = moments_out.avg.avg_raw_s_exiters/moments_out.avg.avg_raw_s;
if moments_out.bhat_probit(3)==0
    xx(5)=1e6; % i.e. as if only s predicts exit
else
    xx(5) = moments_out.bhat_probit(3) / moments_out.bhat_probit(2);
end
xx(6) = moments_out.stats.sd_ik;
xx(7) = moments_out.stats.rho_log_tfpr;
xx(8) = moments_out.stats.sd_log_tfpr;
xx(9) = 1-moments_out.stats.ik_inaction_050; % Pr(|i/k|>x)
wgt = wgt/sum(wgt);

% loss function
qloss = 100*(xx(:)./targets(:) - 1);
qloss = qloss(:).^(2*2);
qloss = qloss'*wgt(:);



