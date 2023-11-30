function output = get_moments_convex(par,N,gr,Px,pol,sim_sw,dat_out_sw)
% Construct all the moments of interest

% keep seed for patternsearch
state_rng = rng;

% set seed
rng('default');

%% Non-simulation Method

% grid for raw mpk
mpk_grid = log(gr.py_i) - log(gr.k_xgrid) + log(par.alpha_k);

% number of quintiles
nq = 3;

% sort MPK
[mpk_raw_s,i_mpk_raw_s] = sort(mpk_grid(:));
% create index for un-sorting
i_inv_sort_mpk_raw(i_mpk_raw_s) = 1:numel(Px.g_ss);

% storage matrix for pr(rw'|rw)
pr_mpk_raw = zeros(nq,nq);

% storage vector for rho_mpk
rho_yklead_yk_vec = zeros(nq,1);

%%% distribution of firms

% initial distribution of active firms
lamb_0 = Px.g_ss(:);
% normalize to measure 1
lamb_0 = lamb_0./sum(lamb_0(:));

% subset dist of continuing firm
lamb_stay = lamb_0(:).*pol.pr_stay(:);
% normalize to 1
lamb_stay = lamb_stay/sum(lamb_stay);

% Next period distribution of firms coming from this period
lamb_1 = Px.H_s*reshape(Px.G_k*lamb_0,N.s,N.k);
% store normalizing measure - needed for normalizing p(x(t)y(t+1)|x(t))*p(x)
m1 = sum(lamb_1(:));
% normalize to 1
lamb_1 = lamb_1(:)./sum(lamb_1(:));

%%%% compute transition matrix

% sort today's distribution (using only stayers)
s_lamb_0_stay = lamb_0(i_mpk_raw_s);
% compute cdf
cs_s_lamb_0_stay = cumsum(s_lamb_0_stay);

% sort tomorrow's distribution
s_lamb_1 = lamb_0(i_mpk_raw_s);
% compute cdf
cs_s_lamb_1 = cumsum(s_lamb_1);

% store average mpk in quantile
mpk_raw_0_q_e = zeros(nq,1);
mpk_raw_1_q_e = zeros(nq,1);
mpk_quant_0 = zeros(nq,1);
mpk_quant_1 = zeros(nq,1);
mpk_grid_q = zeros(nq,length(lamb_0));
E_grk_q = zeros(nq,1);
E_ik_q = zeros(nq,1);
% construct i/k and winsorize
ik_gr = pol.kpr./gr.k_xgrid - (1-par.delta);
    % winsorize window
winsor_w = 0.025;
    % sort ik vector for winsorization
ik_vec = ik_gr(:);
[s_ik_vec,i_s_ik_vec] = sort(ik_vec);
    % sort distribution for winsorization
s_lamb_stay = lamb_stay(i_s_ik_vec);
cs_lamb_stay = cumsum(s_lamb_stay);
    % find window
l_w = find(cs_lamb_stay<=winsor_w,1,'last');
r_w = find(cs_lamb_stay>=(1-winsor_w),1,'first');
    % replace tail values
ik_gr(ik_gr<s_ik_vec(l_w)) = s_ik_vec(l_w);
ik_gr(ik_gr>s_ik_vec(r_w)) = s_ik_vec(r_w);
clear l_w r_w s_ik_vec i_s_ik_vec s_lamb_stay cs_lamb_stay
for iq = 1:nq

    % locate which part of log_yk today corresponds to q (using full dist of
    % today)
    lq = (iq-1)/nq;
    uq = (iq/nq)*(iq<nq) + 100*(iq==nq);
    i_q = cs_s_lamb_0_stay>=lq & cs_s_lamb_0_stay<uq;
    i_q_unsort = i_q(i_inv_sort_mpk_raw);

    % store average r in this quantile
    temp_r = mpk_raw_s(i_q);
    temp_lamb = s_lamb_0_stay(i_q);
    mpk_raw_0_q_e(iq) = temp_r(:)'*temp_lamb(:);
    mpk_quant_0(iq) = find(i_q==1,1,'first');
    clear temp*

    % keep that part of the state vector (original grid)
    lamb_0_i = i_q_unsort.*lamb_0;

    % push one step forward
    lamb_0_i_next = Px.H_s*reshape(Px.G_k*lamb_0_i,N.s,N.k);
    for iq_next = 1:nq
        
        % locate which part of log_yk tomorrow corresponds to quantile q
        lq = (iq_next-1)/nq;
        uq = (iq_next/nq)*(iq_next<nq) + 100*(iq_next==nq);
        i_q_next = cs_s_lamb_1>=lq & cs_s_lamb_1<uq;
        i_q_next_unsort = i_q_next(i_inv_sort_mpk_raw);

        % store average r in this quantile
        temp_r = mpk_raw_s(i_q_next);
        temp_lamb = s_lamb_0_stay(i_q_next);
        mpk_raw_1_q_e(iq_next) = temp_r(:)'*temp_lamb(:);
        mpk_quant_1(iq_next) = find(i_q_next==1,1,'first');
        clear temp*
    
        % compute pr(iq(t+1) | iq(t))
        pr_mpk_raw(iq,iq_next) = lamb_0_i_next(:)'*i_q_next_unsort(:);

    end
    
    % store location on grid for computing regression later
    mpk_grid_q(iq,:) = i_q_unsort(:).*mpk_grid(:);    
    
    % store investment rates in this quantile
    lamb_stay_q = i_q_unsort(:).*lamb_stay(:);
    E_grk_q(iq) = lamb_stay_q(:)'*pol.grk(:)/sum(lamb_stay_q);
    E_ik_q(iq) = lamb_stay_q(:)'*ik_gr(:)/sum(lamb_stay_q);
    
end

output.stats.pr_mpk_raw = pr_mpk_raw;
output.stats.pr_mpk_raw_norm = bsxfun(@rdivide,pr_mpk_raw,sum(pr_mpk_raw,2));
output.stats.rho_yklead_yk = rho_yklead_yk_vec;
output.stats.E_grk_q = E_grk_q;
output.stats.E_ik_q = E_ik_q;

%%% unconditional autocorrelation
% E(yk_(t+1))
E_yklead_all = mpk_grid(:)'*lamb_1(:);
% E(yk(t) | stay)
E_yk_all = mpk_grid(:)'*lamb_stay(:);
% cov(arpk_(t),arpk(t+1))
temp_all = Px.H_s*reshape(Px.G_k*(mpk_grid(:).*lamb_0(:)),N.s,N.k)/m1;
cov_yklead_yk_all = mpk_grid(:)'*temp_all(:) - E_yklead_all*E_yk_all;
% sigma(yk(t+1))
sig_yklead_all = sqrt(sum(((mpk_grid(:)-E_yklead_all).^2).*lamb_1(:)));
% sigma(yk(t) | stay)
sig_yk_all = sqrt(sum(((mpk_grid(:)-E_yk_all).^2).*lamb_stay(:)));
% rho(yk(t+1),yk(t))
output.stats.rho_yklead_yk_all = cov_yklead_yk_all/(sig_yklead_all*sig_yk_all);


%%% Regressions: prob of staying on tfp and k
% turn grids in to logs
log_tfp_grid = log(gr.s_xgrid(:));
log_k_grid = log(gr.k_xgrid(:));
stay_y = pol.pr_stay(:); % split each (s,k) into stay and exit
% use lamb_0 as measure since we use whole popln of active firms
output.stats.beta_hat = pop_reg(stay_y,[log_tfp_grid(:)],lamb_0(:)); % Regression of exit on TFP
output.stats.beta_hat_k =pop_reg(stay_y,[log_k_grid(:)],lamb_0(:)); % Regression of exit on capital
output.stats.beta_hat_all = pop_reg(stay_y,[log_tfp_grid(:) log_k_grid(:)],lamb_0(:)); % Regression of exit on TFP and capital

%%% Regression of mpk(t+1) on mpk(t) in different ranks
% use distribution of continuers for "today's" distribution
output.stats.rho_hat_all = pop_reg(mpk_grid(:),mpk_grid(:),lamb_stay(:),Px,N);
output.stats.rho_hat_condn = pop_reg(mpk_grid(:),mpk_grid_q,lamb_stay(:),Px,N);

%%% compute other moments: growth rate for +ve / -ve shocks
try
    grk_vec = reshape(grk_vec,N.s,N.k);
    grk_vec_low = squeeze(grk_vec(:,1,:));
    grk_vec_high = squeeze(grk_vec(:,2,:));
    lamb_0 = reshape(lamb_0,N.s,N.fstay,N.k);
    lamb_0_low = squeeze(lamb_0(:,1,:));
    lamb_0_high = squeeze(lamb_0(:,2,:));

    output.avg.grk_low_high = [grk_vec_low(:)'*lamb_0_low(:)/sum(lamb_0_low(:)) ; grk_vec_high(:)'*lamb_0_high(:)/sum(lamb_0_high(:))];
catch
    if par.cali == 0
%         disp('only 1 sector')
    end
end

% autocorrelation of TFPR
rhohat = pop_reg(log(gr.s_xgrid(:)),log(gr.s_xgrid(:)),lamb_stay(:),Px,N);
output.stats.rho_log_tfpr = rhohat(2);
clear rhohat

% standard dev of TFPR
lamb_0_new = lamb_0/sum(lamb_0(:));
mean_log_tfpr = log(gr.s_xgrid0(:))'*lamb_0_new(:);
output.stats.sd_log_tfpr = sqrt((log(gr.s_xgrid0(:))'-mean_log_tfpr).^2 * lamb_0_new(:));

% standard dev of MPK
mean_mpk = mpk_grid(:)'*lamb_0(:);
output.stats.sd_mpk = sqrt((mpk_grid(:)'-mean_mpk).^2 * lamb_0(:));

% investment vs MPK
output.stats.beta_hat_grk_mpk = pop_reg(pol.grk(:),[mpk_grid(:)],lamb_stay(:));
output.stats.beta_hat_ik_mpk = pop_reg(ik_gr(:),[mpk_grid(:)],lamb_stay(:));


%% Compute aggregates

%%% Prices
% price of consumption bundle
P = par.P;
Q = par.Q;

%%% Measure of firms
g0  = squeeze(Px.g_ss);                                 % distribution of active firms
g0_entrants = pol.pr_fstar_enter(:).*gr.ps_dist(:);     % distribution of entrants
g0_stay = pol.pr_stay(:).*g0(:);                        % distribution stayers
g0_exit = (1-pol.pr_stay(:)).*g0(:);                    % distribution exiters

%%% Grids
incx_gr = pol.kpr - (1-par.delta)*gr.k_xgrid;    % gross investment grid
incx_gr_true = Q*incx_gr.*(incx_gr>0) + ...
    Q*(1-par.lam)*incx_gr.*(incx_gr<0);          % gross investment grid, factoring in q and Q
ik_i = incx_gr(:)./gr.k_xgrid(:);                   % firm level i/k
kpr_entrants = pol.kpr_up;                          % investment grid of entrants

%%% Aggregates
C = 1/(par.P*par.chi);    % total consumption of final good
% inner integral for domestic goods
inner_Cd = gr.y_i(:).^((par.epsi-1)/par.epsi)'*g0(:);
Cd = inner_Cd^(par.epsi/(par.epsi-1)); % total consumption of domestic good
try 
    Cf = pol.cf_i*par.dM;               % total physical consumption, foreign goods
    P_Cf = par.Pf*pol.cf_i*par.dM;          % total real revenue of foreign firms
catch
    Cf = 0;
    P_Cf = 0;
end
P_Cd = gr.py_i(:)'*g0(:);   % total real revenue of domestic firms
Kd = gr.k_xgrid(:)'*g0(:);  % total capital stock of domestic firms
Ld = gr.lstar(:)'*g0(:);    % total labor of domestic firms
TFPR = P_Cd/(Kd^par.alph * Ld^(1-par.alph)); % TFPR of domestic firms

%%% Investment aggregates (need to account for change in measure)
window = 0.00;
% aggregates
IK = ik_i'*g0_stay(:);                                  % Average I/K of stayers
I_Gross = incx_gr(:)'*g0_stay(:) + ...
        (1-par.delta)*(-1*gr.k_xgrid(:))'*g0_exit + ...
        kpr_entrants(:)'*g0_entrants(:);                % Gross inv incl stayers, exit and entry
I_Gross_plus = (incx_gr(:).*(incx_gr(:)>window))'*g0_stay(:) + ...
        kpr_entrants(:)'*g0_entrants(:);                % Gross +ve incl stayers and entry;
I_Gross_neg = (abs(incx_gr(:)).*(incx_gr(:)<-window))'*g0_stay(:) ...
        + (1-par.delta)*gr.k_xgrid(:)'*g0_exit;      % Gross -ve incl stayers and exits
I_net = incx_gr_true(:)'*g0_stay(:) + ...
        (1-par.zet)*par.Q*(1-par.lam)*(1-par.delta)*(-1*gr.k_xgrid(:)'*g0_exit) + ...
        par.Q*kpr_entrants(:)'*g0_entrants(:);                % Net inv incl stayers, exit, and entry
% Expected values (measures normalized to 1)
I_Gross_intensive = incx_gr(:)'*g0_stay(:)/sum(g0_stay); % Gross inv along intensive margin
I_Gross_plus_intensive = (incx_gr(:).*(incx_gr(:)>window))'*g0_stay(:)/sum(g0_stay); % Gross inv along intensive margin
I_Gross_neg_intensive = (abs(incx_gr(:)).*(incx_gr(:)<-window))'*g0_stay(:)/sum(g0_stay); % Gross inv along intensive margin
Kd_stay = gr.k_xgrid(:)'*g0_stay(:);  % total capital stock of continuing firms
I_Gross_frac_intensive_neg = (incx_gr(:) < -window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of negative investment of stayers
I_Gross_frac_0_intensive = (incx_gr(:)>=-window & incx_gr(:)<=window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of inaction of stayers
I_Gross_frac_pos_intensive = (incx_gr(:) > window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of pos investment of stayers
I_Gross_frac_neg_intensive = (incx_gr(:) < -window)'*g0_stay(:)/sum(g0_stay(:)); % Total fraction of neg investment of stayers

%%% Other costs, labor and welfare
operating_cost = pol.Ef_pr(:)'*g0_stay(:); % total operating costs paid
entry_cost = pol.Ef_enter(:)'*g0_entrants(:); % total entry costs paid
L = Ld + I_net + operating_cost + entry_cost + P_Cf; % total labor (domestic + investment goods + fixed costs + export)
Uss = 1/(1-par.betaa)*(log(C) - par.chi*(L)); % welfare

%%% Raw average of TFPQ
avg_raw_s = gr.s_xgrid0(:)'*g0(:)/sum(g0(:));
avg_raw_s_entrants = gr.s_grid(:)'*g0_entrants(:)/sum(g0_entrants(:));
avg_raw_s_exiters = gr.s_xgrid0(:)'*g0_exit(:)/sum(g0_exit(:));
avg_raw_s_cont = gr.s_xgrid0(:)'*g0_stay(:)/sum(g0_stay(:));

% Measured corrected TFP
M_norm = sum(g0(:)).^(1/(par.epsi-1));
TFPQ = Cd./(M_norm.*(Kd.^par.alph .* Ld.^(1-par.alph)));
    
%%% weighted average of TFPQ (using PY)
PY_domestic = gr.py_i(:)'*g0(:); % PY as denom
wgt_py = (gr.py_i(:).*g0(:))/PY_domestic; % wgts, PY
avg_wgt_s = gr.s_xgrid0(:)'*wgt_py;

%%% weighted average of TFPQ (using Y)
Y_domestic = gr.y_i(:)'*g0(:); % Y as denom
wgt_y = (gr.y_i(:).*g0(:))/Y_domestic; % wgts, Y
avg_wgt_y_s = gr.s_xgrid0(:)'*wgt_y;    

%%% weighted average of TFPQ (K)
K_domestic = gr.k_xgrid(:)'*g0(:); % K as denom
wgt_k = (gr.k_xgrid(:).*g0(:))/K_domestic; % wgts, K
avg_wgt_k_s = gr.s_xgrid0(:)'*wgt_k;   

%%% s.d. TFPQ of continuers
g0_stay_1d = sum(reshape(g0_stay,size(gr.s_xgrid0)),2);
g0_next = g0_stay_1d + g0_entrants;
avg_s = gr.s_grid(:)'*g0_next(:)/sum(g0_next(:));
sd_s_next = sqrt((gr.s_grid(:)'-avg_s).^2 * g0_next(:)/sum(g0_next(:)));

%%% median TFP
g0_s = sum(reshape(g0,N.s,N.k),2);
cdf_g0_s = cumsum(g0_s)/sum(g0_s);
is = find(cdf_g0_s>=0.5,1,'first');
s_med = gr.s_grid(is);
clear g0_s cdf_g0_s

%%% check how many firms exist in "weird" region
% break ties
temps = linspace(sqrt(eps),2*sqrt(eps),N.s);
temps = repmat(temps(:),1,N.k);
pr_stay_new = pol.pr_stay + temps;
[~,ix] = max(pr_stay_new,[],1);
pr0 = 0;
for ik=1:N.k
    g_s = Px.g_ss(ix(ik):end,ik);
    pr0 = sum(g_s(:)) + pr0;
end
frac_odd_region = pr0/sum(Px.g_ss(:));

%%% IK moments
% inaction windows
window_005 = 0.05;
window_010 = 0.10;
window_020 = 0.20;
window_050 = 0.49; % 0.2*sd(LMT)/sd(CH) (sd(i/k))
% normalize distribution of stayers (only stayers in data)
g0_norm = g0_stay(:)./sum(g0_stay(:));
% winsorize "data"
    % winsorize window
winsor_w = 0.025;
    % sort ik vector for winsorization
ik_vec = ik_i(:);
[s_ik_vec,i_s_ik_vec] = sort(ik_vec);
    % sort distribution for winsorization
s_g0_norm = g0_norm(i_s_ik_vec);
cs_g0_norm = cumsum(s_g0_norm);
    % find window
l_w = find(cs_g0_norm<=winsor_w,1,'last');
r_w = find(cs_g0_norm>=(1-winsor_w),1,'first');
    % replace tail values
s_ik_vec_w = s_ik_vec(:);
s_ik_vec_w(1:l_w) = s_ik_vec(l_w);
s_ik_vec_w(r_w:end) = s_ik_vec(r_w);
% compute moments
    % Mean
E_IK = s_ik_vec_w'*s_g0_norm(:);
    % Median
i_med_ik = find(cs_g0_norm>=0.5,1,'first');
med_IK = s_ik_vec_w(i_med_ik);
    % Std Dev
sd_IK = (s_ik_vec_w - E_IK).^2'*s_g0_norm(:);
    % fraction < 0
frac_ik_less_0 = (s_ik_vec_w<0)'*s_g0_norm(:);
    % fraction < -20%
frac_ik_less_20 = (s_ik_vec_w<-window_020)'*s_g0_norm(:);
    % fraction < -29%
frac_ik_less_49 = (s_ik_vec_w<-window_050)'*s_g0_norm(:);
    % E[i/k | ik<0]
E_ik_less_0 = sum(s_ik_vec_w(s_ik_vec_w<0).*s_g0_norm(s_ik_vec_w<0))./frac_ik_less_0;
% E[i/k | ik>0]
frac_ik_pos_000 = (s_ik_vec_w>0.001)'*s_g0_norm(:);
E_ik_pos_000 = sum(s_ik_vec_w(s_ik_vec_w>0.001).*s_g0_norm(s_ik_vec_w>0.001))./frac_ik_pos_000;
% E[i/k | ik>0.20]
frac_ik_pos_020 = (s_ik_vec_w>0.20)'*s_g0_norm(:);
E_ik_pos_020 = sum(s_ik_vec_w(s_ik_vec_w>0.20).*s_g0_norm(s_ik_vec_w>0.20))./frac_ik_pos_020;
% E[i/k | ik>0.5]
frac_ik_pos_050 = (s_ik_vec_w>0.5)'*s_g0_norm(:);
E_ik_pos_050 = sum(s_ik_vec_w(s_ik_vec_w>0.5).*s_g0_norm(s_ik_vec_w>0.5))./frac_ik_pos_050;
% E[i/k | ik>0.75]
frac_ik_pos_075 = (s_ik_vec_w>0.75)'*s_g0_norm(:);
E_ik_pos_075 = sum(s_ik_vec_w(s_ik_vec_w>0.75).*s_g0_norm(s_ik_vec_w>0.75))./frac_ik_pos_075;
    % inaction regions
ik_inaction_005 = (abs(s_ik_vec_w)<=window_005)'*s_g0_norm;
ik_inaction_010 = (abs(s_ik_vec_w)<=window_010)'*s_g0_norm;
ik_inaction_020 = (abs(s_ik_vec_w)<=window_020)'*s_g0_norm;
ik_inaction_050 = (abs(s_ik_vec_w)<=window_050)'*s_g0_norm;
% compile into a table (as in appendix table A1)
table_ik_moments = [E_IK med_IK sd_IK frac_ik_less_0 E_ik_less_0 ik_inaction_005 ik_inaction_010 ik_inaction_020 ik_inaction_050 E_ik_pos_000 E_ik_pos_020 E_ik_pos_050 E_ik_pos_075];
output.stats.sd_ik = sd_IK;
output.stats.ik_inaction_050 = ik_inaction_050;
output.stats.frac_ik_less_20 = frac_ik_less_20;
output.stats.frac_ik_less_49 = frac_ik_less_49;

%%% Compact into single output
output.avg.P = P;
output.avg.Q = Q;
output.avg.C = C;
output.avg.Cd = Cd;
output.avg.Cf = Cf;
output.avg.P_Cd = P_Cd;
output.avg.P_Cf = P_Cf;
output.avg.Kd = Kd;
output.avg.Kd_stay = Kd_stay;
output.avg.Ld = Ld;
output.avg.L = L;
output.avg.TFPR = TFPR;
output.avg.TFPQ = TFPQ;
output.avg.Mentrants = sum(g0_entrants(:));
output.avg.Mactive = sum(g0(:));
output.avg.Mexit = sum(g0_exit(:));
output.avg.exit_rate = output.avg.Mexit/output.avg.Mactive;

output.avg.IK = IK;
output.avg.I_Gross = I_Gross;
output.avg.I_Gross_plus = I_Gross_plus;
output.avg.I_Gross_neg = I_Gross_neg;
output.avg.I_Gross_frac_0_intensive = I_Gross_frac_0_intensive;
output.avg.I_Gross_frac_pos_intensive = I_Gross_frac_pos_intensive;
output.avg.I_Gross_frac_neg_intensive = I_Gross_frac_neg_intensive;
output.avg.I_net = I_net;
output.avg.I_Gross_intensive = I_Gross_intensive;
output.avg.I_Gross_plus_intensive = I_Gross_plus_intensive;
output.avg.I_Gross_neg_intensive = I_Gross_neg_intensive;
output.avg.I_Gross_frac_neg = I_Gross_frac_intensive_neg;

output.avg.operating_cost = operating_cost;
output.avg.entry_cost = entry_cost;
output.avg.Uss = Uss;

output.avg.avg_raw_s = avg_raw_s;
output.avg.avg_wgt_s = avg_wgt_s;
output.avg.avg_wgt_y_s = avg_wgt_y_s;
output.avg.avg_wgt_k_s = avg_wgt_k_s;

output.avg.avg_raw_s_entrants = avg_raw_s_entrants;
output.avg.avg_raw_s_exiters = avg_raw_s_exiters;
output.avg.avg_raw_s_cont = avg_raw_s_cont;

output.avg.s_med = s_med;

output.avg.sd_s_next = sd_s_next;

output.avg.frac_odd_region = frac_odd_region;

output.avg.table_ik_moments = table_ik_moments;

%%% pull out moments
% average k of current active firms
output.moments.avg_k_all = gr.k_xgrid(:)'*g0(:)/sum(g0(:));
% average k of exiters
output.moments.avg_k_exit = gr.k_xgrid(:)'*g0_exit(:)/sum(g0_exit(:));
% average k of continuers
output.moments.avg_k_cont = gr.k_xgrid(:)'*g0_stay(:)/sum(g0_stay(:));
% average k of entrants
output.moments.avg_k_enter = kpr_entrants(:)'*g0_entrants(:)/sum(g0_entrants(:));
% probs k'/k<0 for continuers
neg_kpr_gr = pol.kpr < gr.k_xgrid;
output.moments.avg_frac_neg_kpr = neg_kpr_gr(:)'*g0_stay(:)/sum(g0_stay(:));
% probs k' < (1-delta)k for continuers
output.moments.I_Gross_frac_intensive_neg = I_Gross_frac_intensive_neg;

%% Simulation
% construct MPK using simulation
if sim_sw==1 && output.avg.exit_rate>0.01 % do only if exit rate >1%
    
    % simulation parameters
    T = 50;
    t_burn = 30;
    Mactive=sum(g0(:));

    %%% Grid policy functions (in linear form)
    % continuation kprime policy
    FK_cont = griddedInterpolant({1:N.s,gr.k_grid},pol.kpr,'linear');
    % exit probs, conditioned on state (s,k)
    Fexit = griddedInterpolant({1:N.s,gr.k_grid},pol.pr_stay,'linear');
    % entry probs, conditioned on state (s)
    Fentry = griddedInterpolant({1:N.s},pol.pr_fstar_enter,'linear');
    % entry kprime policy
    FK_entry = griddedInterpolant({1:N.s},pol.kpr_up,'linear');
    % exit choice, conditioned on state (s,k)
    Fstar_exit = griddedInterpolant({1:N.s,gr.k_grid},pol.fstar,'linear');
    % entry choice, conditioned on state (s)
    Fstar_enter = griddedInterpolant({1:N.s},pol.f_star_enter,'linear');
    
    % actual grid numbers
    [i_s_xgrid , k_xgrid] = ndgrid(1:N.s,gr.k_grid);
    i_s_xgrid = i_s_xgrid(:);
    k_xgrid = k_xgrid(:);
    
    %%% construct initial panel of active firms
    % frequency at each state (s,fstay,k)
    prk = round(Px.g_ss(:)/sum(Px.g_ss(:))*10000);
    c_prk = cumsum(prk);
    N_active_star = sum(prk);
    % total length of panel
    N_active = sum(prk); 
    % store initial conditions
    k_active = zeros(N_active,1);
    i_s_active = zeros(N_active,1);
    for ii=1:length(prk)
        if ii==1
            istart = 1;
        else
            istart = c_prk(ii-1)+1;
        end
        iend = c_prk(ii);
        
        k_active(istart:iend) = repmat(k_xgrid(ii),prk(ii),1);
        i_s_active(istart:iend) = repmat(i_s_xgrid(ii),prk(ii),1);
    end
    % invariant distribution of s
    ps_pdf = gr.Ps^1000;
    ps_pdf = ps_pdf(1,:);

    
    kcont_sim = cell(T,1);
    kenter_sim = cell(T,1);
    i_scont_sim = cell(T,1);
    i_senter_sim = cell(T,1);
    N_active_sim = cell(T,1);
    i_enter_sim = cell(T,1);
    i_stay_sim = cell(T,1);
    k_active_sim = cell(T,1);
    i_s_active_sim = cell(T,1);
    i_track_sim = cell(T,1);
    i_s_active_enter_sim = cell(T,1);
    
    % save all initial conditions
    kcont_sim{1} = k_active;
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
    G_f_s = @(fm,sgrid) fm.*(par.s_sh - (1+sgrid).^par.eta(1) + (1+sgrid).^(par.eta(1)+par.eta(2)));
    % Inverse cdf for sampling f
    igc = @(pr,cm) pr.^(1/par.npow).*cm;
    %%%
    
%     %%% Version 1: Original with mixture probs
%     % Simulate using ver 0 functions, but set fstar=Inf for x% of the time
%     % (i.e. for x% of the time, firms will not exit
%     % ex: f_pot_stay = igc(rand(N_active,1),xm).*(rand(N_active,1)<par.s_sh) + 1e10.*(rand(N_active,1)>=par.s_sh);
%     %%%
    
%     %%% Version 2: g(f,fstar)
%     % the probability function is a function of fstar itself. This allows
%     % direct control of the slope of Pr(stay,s)
%     % G(f,s) function: Does not directly depend on f or s
%     G_f_s = @(fm,sgrid) fm;
%     % Inverse cdf for sampling f
%     igc = @(pr,cm) (pr*(cm/sx_sh)).^(1/(k-1/sk(2)));
%     %%%
    
    for t=2:T
    
        %%%%%%%%% Current cohort's choices %%%%%%%%%

        N_active = N_active_sim{t-1};
        i_s_active = [i_scont_sim{t-1}(:) ; i_senter_sim{t-1}(:)];
        k_active = [kcont_sim{t-1}(:) ; kenter_sim{t-1}(:)];
        i_track = [0*kcont_sim{t-1}(:) ; 1+0*kenter_sim{t-1}(:)];

        % draw current period fstay 
        xm = G_f_s(par.mu_fstay_b,gr.s_grid(i_s_active).^((par.epsi-1)/par.epsi));
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
        N_pot_entrants = round((par.Me/Mactive)*N_active_star);

        % draw current period s for potential entrants (iid!)
        i_s_pot_enter = sim_discrete_iid(ps_pdf,N_pot_entrants,'pdf');
        i_s_pot_enter = i_s_pot_enter(:);
        
        % draw current period fstay 
        xm = G_f_s(par.mu_fenter_b,gr.s_grid(i_s_pot_enter).^((par.epsi-1)/par.epsi));
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

        %%% save all variables
        
        % choices 
        kcont_sim{t} = kprime_cont;
        kenter_sim{t} = kprime_enter;
        i_scont_sim{t} = i_sprime_stay;
        i_senter_sim{t} = i_sprime_enter;
        N_active_sim{t} = Nstay+Nenter;
        i_enter_sim{t} = i_enter;
        i_stay_sim{t} = i_stay;
        
        % states
        k_active_sim{t} = k_active;
        i_s_active_sim{t} = i_s_active;
        i_track_sim{t} = i_track;
        i_s_active_enter_sim{t} = [i_s_active(i_stay) ; i_s_enter];
        
    end
    
    %%% pull out data 
    % drop burn in data
    kcont_sim = kcont_sim(t_burn:end-1);
    i_scont_sim = i_scont_sim(t_burn:end-1);
    k_active_sim = k_active_sim(t_burn:end);
    i_s_active_sim = i_s_active_sim(t_burn:end);
    i_s_active_enter_sim = i_s_active_enter_sim(t_burn-1:end-1);
    i_stay_sim = i_stay_sim(t_burn:end);
    % convert to single vector
    kpr_vec = cell2mat(kcont_sim);
    i_spr_vec = cell2mat(i_scont_sim);
    k_vec = cell2mat(k_active_sim);
    i_s_vec = cell2mat(i_s_active_sim);
    i_stay_vec = cell2mat(i_stay_sim);
    T2 = T - t_burn;
    % pull out actual k and s values
    k_sim = k_vec;
    s_sim = gr.s_grid(i_s_vec).^((par.epsi-1)/par.epsi);
    
    if sum(g0_exit(:))/sum(g0(:))>.001
        %%% linear regression
        output.bhat_ols = regress(i_stay_vec,[ones(size(i_s_vec)) log(s_sim) log(k_sim)]);

        %%% probit
        output.bhat_probit = glmfit([log(s_sim) log(k_sim)],i_stay_vec,'binomial','link','probit');
    else
        %%% linear regression
        output.bhat_ols = regress(i_stay_vec,[ones(size(i_s_vec)) log(s_sim) log(k_sim)]);

        %%% probit
        output.bhat_probit = [0 0 0];
    end
    
    %%% transition matrix
    nq = 3;
    logmpk_stay_sim_qvec = zeros(nq+1,T2-1);
    logmpk_stay_prime_sim_qvec = zeros(nq+1,T2-1);
    logmpk_sim_qmat = cell(T2-1,1);
    logmpk_prime_sim_qmat = cell(T2-1,1);
    mpk_now = cell(T2-1,1);
    mpk_tmr = cell(T2-1,1);
    for t=1:T2-1
        
        % compute MPK for every active firm
        logmpk_sim = log(par.alpha_k) + log(gr.s_grid(i_s_active_sim{t}(:)).^((par.epsi-1)/par.epsi)) + (par.alpha_k-1)*log(k_active_sim{t}(:));
        
        % subset out firms that are continuing
        logmpk_stay_sim = logmpk_sim(i_stay_sim{t});
        mpk_now{t} = logmpk_stay_sim;
        
        % construct quantiles for these guys
        n_bins = 2;
        noise_correction = 1e-10;
        logmpk_stay_sim2 = logmpk_stay_sim + (rand(size(logmpk_stay_sim)) - 0.5).*noise_correction;
        logmpk_stay_sim_qvec(:,t) = [min(logmpk_stay_sim2) quantile(logmpk_stay_sim2',n_bins) max(logmpk_stay_sim2)];
        logmpk_sim_qmat{t} = discretize(logmpk_stay_sim2,logmpk_stay_sim_qvec(:,t));
        
        %%% Quantiles for tomorrow
        
        % construct quantiles for continuers into tomorrow
        kcont = kcont_sim{t}; %(1:sum(i_stay_sim{t}));
        scont = gr.s_grid(i_scont_sim{t}).^((par.epsi-1)/par.epsi); %(1:sum(i_stay_sim{t}));
        
        % compute MPK for every active firm
        logmpk_prime_sim = log(par.alpha_k) + log(scont(:)) + (par.alpha_k-1)*log(kcont(:));
        mpk_tmr{t} = logmpk_prime_sim;

        % construct quantiles for these guys
        n_bins = 2;
        noise_correction = sqrt(eps);
        logmpk_prime_sim2 = logmpk_prime_sim + (rand(size(logmpk_prime_sim)) - 0.5).*noise_correction;
        try
            switch q_set
                case('t,t+1')
                    logmpk_stay_prime_sim_qvec(:,t) = [min(logmpk_prime_sim2) quantile(logmpk_prime_sim2',n_bins) max(logmpk_prime_sim2)];
                    logmpk_prime_sim_qmat{t} = discretize(logmpk_prime_sim2,logmpk_stay_prime_sim_qvec(:,t));
                case('t,t')
                    logmpk_stay_prime_sim_qvec(:,t) = [min(logmpk_prime_sim2) logmpk_stay_sim_qvec(2:3,t)' max(logmpk_prime_sim2)];
                    logmpk_prime_sim_qmat{t} = discretize(logmpk_prime_sim2,logmpk_stay_prime_sim_qvec(:,t));
            end
        catch
            % default, use a single fixed quantile
            logmpk_stay_prime_sim_qvec(:,t) = [min(logmpk_prime_sim2) logmpk_stay_sim_qvec(2:3,t)' max(logmpk_prime_sim2)];
            logmpk_prime_sim_qmat{t} = discretize(logmpk_prime_sim2,logmpk_stay_prime_sim_qvec(:,t));
        end
        
    end
    % bin and compute transition probabilities
    trq = zeros(T2-1,nq,nq);
    for t = 1:T2-1
        one = logmpk_sim_qmat{t};
        two = logmpk_prime_sim_qmat{t};
        for qq1 = 1:nq
            for qq2 = 1:nq
                trq(t,qq1,qq2) = sum(one == qq1 & two == qq2)/sum(one==qq1 & isfinite(two));
            end
        end
    end
    output.pr_logmpk_hat = squeeze(mean(trq,1));

    %%% AR(1) regressions
    i_mpk_now = cell2mat(logmpk_sim_qmat);
    mpk_now = cell2mat(mpk_now);
    mpk_tmr = cell2mat(mpk_tmr);
    X = [ones(size(mpk_now(:))) mpk_now(:).*(i_mpk_now(:)==1) mpk_now(:).*(i_mpk_now(:)==2) mpk_now(:).*(i_mpk_now(:)==3)];
    X_all = [ones(size(mpk_now(:))) mpk_now(:)];
    Y = mpk_tmr(:);
    % AR(1)
    output.bhat_all = regress(Y,X_all);
    % with dummies
    output.bhat = regress(Y,X);
    yhat1 = X*output.bhat;
    
    % investment elasticity regresions
    output.ik_regs = ik_elas(nq,i_mpk_now,mpk_now,kpr_vec,k_active_sim,i_s_active_sim,i_s_active_enter_sim,i_stay_sim,gr,par,T2);
    
    if nargin==7
        if dat_out_sw==1
            output.dat.mpk_now = mpk_now;
            output.dat.mpk_tmr = mpk_tmr;
            output.dat.i_mpk_now = i_mpk_now;
            output.dat.logmpk_stay_prime_sim_qvec = logmpk_stay_prime_sim_qvec;
            output.dat.yhat1 = yhat1;
        end
    end
    
else
    % if exit rate is very small, just straight to here
    
    output.bhat_ols = [0 0 0];
    output.bhat_probit = [0 0 0];
    
end

% restore state of rand
rng(state_rng);

end

function outregs = ik_elas(nq,i_mpk_now,mpk_now,kpr_vec,k_active_sim,i_s_active_sim,i_s_active_enter_sim,i_stay_sim,gr,par,T2)

    % choice of winsorization
    q_winsor = .99;
    
    %%%%%%%%%%%%%%%% Set up simulated data %%%%%%%%%%%%%%%%

    % drop last period (did not compute i_mpk_now)
    kpr_vec_2 = kpr_vec(1:length(i_mpk_now));

    % pull out current k for same period
    k_active_vec = [];
    for tt = 1:T2-1
        k_active_temp = k_active_sim{tt}(i_stay_sim{tt});
        k_active_vec = [k_active_vec ; k_active_temp(:)];
    end
    % log k
    logk = log(k_active_vec);
        % compute growth rate of capital, winsorized
    log_grk =  log(kpr_vec_2./k_active_vec);
    log_grk(log_grk>quantile(log_grk,q_winsor)) = quantile(log_grk,q_winsor);
    log_grk(log_grk<quantile(log_grk,1-q_winsor)) = quantile(log_grk,1-q_winsor);
        % compute investment rate, winsorized
    ik = kpr_vec_2./k_active_vec - (1-par.delta);
    ik(ik>quantile(ik,q_winsor)) = quantile(ik,q_winsor);
    ik(ik<quantile(ik,1-q_winsor)) = quantile(ik,1-q_winsor);

    % pull out lag and current s for same period
    i_s_now_vec = [];
    i_s_lag_vec = [];
    for tt = 1:T2-1
        i_s_active_temp = i_s_active_sim{tt}(i_stay_sim{tt});
        i_s_active_enter_temp = i_s_active_enter_sim{tt}(i_stay_sim{tt});

        i_s_now_vec = [i_s_now_vec ; i_s_active_temp(:)];
        i_s_lag_vec = [i_s_lag_vec ; i_s_active_enter_temp(:)];
    end
    s_now_sim = gr.s_grid(i_s_now_vec).^((par.epsi-1)/par.epsi);
    s_lag_sim = gr.s_grid(i_s_lag_vec).^((par.epsi-1)/par.epsi);
        % compute TFP "shock" as AR(1) residual
    tfp_shock = log(s_now_sim) - par.rho_s*log(s_lag_sim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Regressions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% set up regressors
    % MPK with het intercept
        % note: "X" in orig. code forces common intercept
    X = zeros(length(mpk_now),nq+1);
    X(:,1) = ones(size(mpk_now));
    for ii=1:nq
        X(:,ii+1) = mpk_now.*(i_mpk_now==ii);
    end
    X2 = zeros(length(mpk_now),nq+nq);
    for ii=1:nq
        if ii==1
            X2(:,ii) = X(:,1);
        else
            X2(:,ii) = (i_mpk_now(:)==ii);
        end
    end
    X2(:,ii+1:end) = X(:,2:end);
    % demeaned MPK
    mpk_now_demean = mpk_now - mean(mpk_now(:));
    X_mpk_demean = zeros(length(mpk_now),nq+1);
    X_mpk_demean(:,1) = ones(size(mpk_now_demean(:)));
    for ii=1:nq
        X_mpk_demean(:,ii+1) = mpk_now_demean(:).*(i_mpk_now(:)==ii);
    end
    % TFP with common intercept
    X_tfp = zeros(length(mpk_now),nq+1);
    X_tfp(:,1) = ones(size(tfp_shock(:)));
    for ii=1:nq
        X_tfp(:,ii+1) = tfp_shock(:).*(i_mpk_now(:)==ii);
    end
    % TFP with het intercept
    X2_tfp = zeros(length(mpk_now),nq+nq);
    for ii=1:nq
        if ii==1
            X2_tfp(:,ii) = X_tfp(:,1);
        else
            X2_tfp(:,ii) = (i_mpk_now(:)==ii);
        end
    end
    X2_tfp(:,ii+1:end) = X_tfp(:,2:end);

    %%% run regressions
    
    % univariate uncondition regression, MPK
    [outregs.bhat_grk_mpk_simple,outregs.ci_grk_mpk_simple]=regress(log_grk(:),[ones(size(mpk_now(:))) mpk_now(:)]);
    % common intercept MPK
    [outregs.bhat_grk_mpk,outregs.ci_grk_mpk] = regress(log_grk(:),X);
    % het intercept MPK
    [outregs.bhat_grk_mpk_het,outregs.ci_grk_mpk_het] = regress(log_grk(:),X2);
    % demeaned MPK, common intercept
    [outregs.bhat_grk_mpk_demean,outregs.ci_grk_mpk_demean] = regress(log_grk(:),X_mpk_demean);
    % univariate uncondition regression, TFP shock
    [outregs.bhat_grk_tfp_simple,outregs.ci_grk_tfp_simple]=regress(log_grk(:),[ones(size(tfp_shock(:))) tfp_shock(:)]);
    % TFP shock, common intercept
    [outregs.bhat_grk_tfp,outregs.ci_grk_tfp] = regress(log_grk(:),X_tfp);
    % TFP shock, het intercept
    [outregs.bhat_grk_tfp_het,outregs.ci_grk_tfp_het] = regress(log_grk(:),X2_tfp);

    % inaction ~ intercept_mpk + beta_mpk*mpk
        % LP, univariate uncondition regression, MPK
    [outregs.bhat_i_inact_mpk_simple,outregs.ci_i_inact_mpk_simple]=regress(~(abs(ik)>sqrt(eps)),[ones(size(mpk_now(:))) abs(mpk_now_demean(:))]);
        % LP, het intercept MPK
    [outregs.bhat_i_inact_mpk,outregs.ci_i_inact_mpk] = regress( ~(abs(ik)>sqrt(eps)), X2);
%         % probit, univariate uncondition regression, MPK
%     outregs.bhat_i_inact_simple_probit = glmfit(mpk_now(:),~(abs(ik)>sqrt(eps)),'binomial','link','probit');
%         % probit, het intercept MPK
%     outregs.bhat_i_inact_probit = glmfit(X2(:,2:end),~(abs(ik)>sqrt(eps)),'binomial','link','probit');


end