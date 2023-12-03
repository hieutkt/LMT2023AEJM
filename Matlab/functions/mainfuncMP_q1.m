function output = mainfuncMP_q1(par,N,grids,pol,sw,initval)
% Solves q=Q only

    if sw==1
        % get baseline parameters
        try
            [output.par,output.N]=get_params(initval);
        catch
            [output.par,output.N]=get_params;
        end
        
    elseif sw==2
        % set the grids
        output=setgrid(par,N);
    
    elseif sw==3
        % solve vfi
        try
            [output.Ev,output.pol]=solve_vfi(par,N,grids,initval);
        catch
            [output.Ev,output.pol]=solve_vfi(par,N,grids);
        end
        
    elseif sw==4
        % get s.s. distribution
        [output.g_ss,output.G_k,output.H_s] = get_dist(N,grids,pol);
        
    elseif sw==5
        % one step backward induction for transition
        [output.vf,output.pol]=solve_vfi_onestep(par,N,grids,initval);
        
    elseif sw==6
        % forward iteration associated with this policy
        [output.g1,output.G_k,output.H_s] = get_trmat(N,grids,pol,initval.trmat,initval.g0);
        
    end

end

%% Fix baseline parameters
function [par,N] = get_params(Nnew)

try
    % resale price of used capital
    par.lam = Nnew.lam;                

    % fixed cost: stay
    par.mu_fstay_a = Nnew.mu_fstay_a;
    par.mu_fstay_b = Nnew.mu_fstay_b;
    
    % measure of entrants
    par.Me = Nnew.Me;
    
    % exit cost
    par.zet = Nnew.zet; % exit cost
    
    if Nnew.cali == 0
        disp('varying parameters')
        disp('')
    end
    
catch
    % no parameters supplied
    
    % resale price of used capital
    par.lam = .74;                

    % fixed cost: stay
    par.mu_fstay_a = 1;
    par.mu_fstay_b = 1+0.9;
    
    % measure of entrants
    par.Me = 1;
    
    % exit cost
    par.zet = 0; % exit cost
    
end

    %%% fixed parameters
    par.chi     = 2.15;         % disutility of leisure
    par.epsi    = 4;            % Elasticity of sub btw types
    par.nu      = (par.epsi-1)/par.epsi; % effective returns to scale
%     par.alpha_k    = 0.307;         % actual k share
    par.alpha_k    = 0.297;         % actual k share
    par.alph    = par.alpha_k/par.nu;  % capital intensity
    par.alpha_l = par.nu-par.alpha_k; % actual labor share
    par.betaa = .96;             % discount factor
%     par.delta = .11;           % depreciation
    par.delta = .105;           % depreciation
    par.f = 0;                  % fixed operating cost
    par.c = 0;                  % utilization cost
%     par.sigma_us = .812;         % conditional std dev of TFP (transform with scaling!)
%     par.rho_s = .729;             % autocorrelation of TFP
    par.sigma_us = .758;         % conditional std dev of TFP (transform with scaling!)
    par.rho_s = .742;             % autocorrelation of TFP
    par.sigma_s = par.sigma_us/sqrt(1-par.rho_s^2);
    par.mu_shift = -par.sigma_s^2/2;
    par.A = 1; % aggregate TFP
    try
        par.mu_fenter_a = Nnew.mu_fenter_a;
        par.mu_fenter_b = Nnew.mu_fenter_b;
    catch
        par.mu_fenter_a = par.mu_fstay_a;
        par.mu_fenter_b = par.mu_fstay_b;
    end
    
    %%% extra new stuff
    % exponent / pareto scale
    par.npow = Nnew.npow;
    % G(f,s) correlation
    par.eta = Nnew.eta;
    % mean shifter
    par.s_sh = Nnew.s_sh;
    % scrappage function elasticity
    par.kapp = Nnew.kapp;
    
    %%% When calibrating sigma(s)
    try
        par.sigma_us = Nnew.sigma_us;
    catch
        disp('Use default sigma')
    end
    
    %%% When calibrating rho(s)
    try
        par.rho_s = Nnew.rho_s;
    catch
        disp('Use default rho')
    end
    
    % num of nodes
    try
        N.s = Nnew.Ns;
    catch
        N.s = 17;
        if Nnew.cali ~= 1
            disp('Use default Ns');
        end
    end
    try
        N.k = Nnew.Nk;
    catch
        N.k = 5000;
        if Nnew.cali ~= 1
            disp('Use default Nk');
        end
    end
    try
        N.fstay = Nnew.fstay;
    catch
        N.fstay = 21;
        if Nnew.cali ~= 1
            disp('Use default Nfstay');
        end
    end
    
    % check if using fixed R or vary R
    try
        par.fixedR = Nnew.fixedR ;
    catch
        % default is GE (vary R)
        par.fixedR = 0;
    end
    
    % assign for later
    par.cali = Nnew.cali;
    
    
    try
        par.C_f = Nnew.C_f;
        par.model_fixed_cost = Nnew.model_fixed_cost;
    catch
        par.C_f = 0;
    end

end

%% Set grids
function grids = setgrid(par,N)

    v2struct(par)
    
    %%% demand externality ("today")
    Yd = 1/(P*chi);
    Afac = A*P*Yd^(1/epsi);

    %%% demand externality ("tomorrow")
        % for MIT shocks
        % Note: Afac_pr = Afac if steady-state!
    try
        Yd_pr = 1/(P_pr*chi);
        Afac_pr = A*P_pr*Yd_pr^(1/epsi);
    catch
        Yd_pr = Yd;
        Afac_pr = Afac;
    end

    %%% s shocks
    [logs_grid,grids.Ps] = tauchen(N.s,mu_shift,rho_s,sigma_us,2.5);
    grids.s_grid = exp(logs_grid)*A;

    %%% invariant distribution of s
    ps_dist = grids.Ps^1000;
    grids.ps_dist = par.Me*ps_dist(1,:);
    
    %%% k grids (no need to set)
    % q=Q means no need to iterate on bellman
    % Note: for ss, k(t+1) and k(t) are the same grid

    % set discounting for analytical solution
    if par.fixedR==1
        % Note: this is never used except for transitional dynamics
        betaa_pr = betaa*P/P_pr;
        betaa_last = betaa*P_last/P;
        Q_last = Q;
    else
        betaa_pr = betaa;
        betaa_last = betaa;
        Q_last = Q;
    end

    % override original grid size, replace directly with "on-grid" size
        % 1 k to 1 s
    N.k = N.s;

    % compute E[s'|s]
    ss = grids.s_grid.^(((epsi-1)/epsi)*(1+alpha_l/(1-alpha_l)));
    Es = grids.Ps*ss;

    % y-wl = f(alp,w)*stuff
    w_diff = alpha_l^(alpha_l/(1-alpha_l))-alpha_l^(1/(1-alpha_l));

    % k(t+1): directly obtain the optimal kpr
    vartheta_pr = Q*(1/betaa_pr-1+delta)*((1-alpha_l)/alpha_k) * Afac_pr.^-(1+alpha_l/(1-alpha_l)) / w_diff;
    grids.kstar_noirr_pr = (vartheta_pr./Es).^(1/(alpha_k-1+alpha_k*alpha_l/(1-alpha_l)));
    % directly back out grid number (trivially 1:Nk)
    grids.i_kpr = 1:N.k;
    grids.i_s_k = sub2ind([N.s,N.k],(1:N.s)',grids.i_kpr(:));

    % kpr(t): Use this to construct the "today" grid
%         vartheta = Q_last*(1/betaa_last-1+delta)/alpha_k * Afac.^-(1+alpha_l/(1-alpha_l)) * alpha_l^-(alpha_l/(1-alpha_l));
    vartheta = Q_last*(1/betaa_last-1+delta)*((1-alpha_l)/alpha_k) * Afac.^-(1+alpha_l/(1-alpha_l)) / w_diff;
    grids.kstar_noirr = (vartheta./Es).^(1/(alpha_k-1+alpha_k*alpha_l/(1-alpha_l)));
    % override original grid, and replace directly with "on-grid" soln
    grids.k_grid = grids.kstar_noirr;
    
    %%% for transitional dynamics, override with ss grid for k at t=0
    try
        if t==1
            grids.k_grid = par.k_grid;
            % disp('t=1')
        end
    catch
%             disp('solve stationary model')
    end

    %%% full mesh grid (needed to compute net y_cont)
    % primary grids
    [grids.s_xgrid0, grids.k_xgrid] =ndgrid(grids.s_grid,grids.k_grid);
    grids.s_xgrid = grids.s_xgrid0.^((epsi-1)/epsi);
    

    % optimal labor
    grids.lstar = (alpha_l*Afac.*grids.s_xgrid.*grids.k_xgrid.^alpha_k).^(1/(1-alpha_l));    
    % optimal profits (before investment)
    grids.y = Afac.*grids.s_xgrid.*grids.k_xgrid.^alpha_k.*grids.lstar.^alpha_l - grids.lstar;
    % nominal revenue
    grids.py_i = Afac.*grids.s_xgrid.*grids.k_xgrid.^alpha_k.*grids.lstar.^alpha_l;
    % physical output
    grids.y_i = grids.s_xgrid0.*grids.k_xgrid.^alph.*grids.lstar.^(1-alph);

    % reshape into exo-endo
    grids.y = reshape(grids.y,N.s,N.k);
    grids.py_i = reshape(grids.py_i,N.s,N.k);
    grids.y_i = reshape(grids.y_i,N.s,N.k);
    grids.k_xgrid = reshape(grids.k_xgrid,N.s,N.k);
    
end

%% VFI solver
function [Evf,pol] = solve_vfi(par,N,grids,initvals)

    % unpack parameters
    v2struct(par)
    v2struct(grids)

    % initial expected value functions
    if nargin==4
        Ev = initvals.Ev;
    else
        Ev = ones(N.s,N.k);
    end    
    % other loop parameters
    diff_v=1;
    maxiter_v = 1000;
    tol_v = 10^-5;
    iter_v = 0;

    % scrappage recovery function
    kapp  = par.kapp + (1-par.kapp)*(zet==0 && lam==0);
    hk = @(ks) (1-zet)*(1-lam)*(1+ks).^(kapp-1);
    

    while diff_v > tol_v && iter_v < maxiter_v

        iter_v = iter_v +1;
        Ev0 = Ev;

        % no resale frictions: directly use analytical solution to compute v
        v_cont = -Q*grids.kstar_noirr + betaa*Ev(i_s_k);

        % pull out value of staying, gross of fixed operating cost
        % (v_stay(s,k))
        v_stay = y + Q*(1-delta).*k_xgrid + repmat(v_cont(:),1,N.k);
        v_up_entry = v_cont; % continuation value for entering
        
        % value of exiting
        v_exit = y + hk(k_xgrid)*Q*(1-delta).*k_xgrid;
        
        % implied f that makes v_stay>v_exit
        f_star = max(sqrt(eps),v_stay - v_exit);
            % Pr(f'<fstar(s',k'))
        [pr_fstar,Ef_pr] =  cdfu(f_star,mu_fstay_a,mu_fstay_b,npow,eta,s_sh,s_xgrid);
            % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
        lhs = (v_stay - Ef_pr).*pr_fstar;
            % (pi(s',k') + Q*(1-lam)*(1-delta)*k')*(1-Pr(f'<fstar(s',k')))
        rhs = v_exit.*(1-pr_fstar);
            % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
        Ev = Ps*(lhs+rhs);

        % check pr>=0
        if min(pr_fstar(:)<0)
            disp('error')
            break
        end

        % compute difference
        diff_v = max(abs(Ev(:)-Ev0(:)));
        
    end
    
    % value of entry
    f_star_enter = v_up_entry;
    [pr_fstar_enter,Ef_enter] =  cdfu(f_star_enter,mu_fenter_a,mu_fenter_b,npow,eta,s_sh,s_grid.^((epsi-1)/epsi));
    
    %%% pull out all policy functions
    pol.i_kpr = repmat(i_kpr(:),1,N.k);
    pol.i_kpr_up = i_kpr;
    pol.kpr = repmat(grids.kstar_noirr(:),1,N.k);
    pol.kpr_up = grids.kstar_noirr;
        % Policy functions common to both lam=0 and lam>0
    pol.grk = log(pol.kpr(:)./k_xgrid(:));
    pol.fstar = f_star;
    pol.f_star_enter = f_star_enter;
    pol.pr_stay = pr_fstar; % transition probs for stayers
    pol.pr_fstar_enter = pr_fstar_enter; % entry prop conditioned on z draw
    pol.Ef_pr = Ef_pr;
    pol.Ef_enter = Ef_enter;
        % value functions
    pol.v_exit = v_exit;
    pol.v_stay = v_stay;
    
    Evf = Ev;    
    
    % check that distribution is not degenerate
    if min(pr_fstar(:)==1)
        disp('warning: exit rate is 0')
        return
    end
    cdf_pr_s = Ps^100;
    cdf_pr_s = cdf_pr_s(1,:);
    if pr_fstar_enter(:)'*cdf_pr_s(:)<=0.01
        disp('warning: entry rate is approximately 0!')
        return
    end
end

%% Solve for s.s. dist
function [g_ss,G_k,H_s] = get_dist(N,gr,pol)

    % total # of states
    n_states = N.k*N.s;
    
    
    
    % I(k',s,fs | k,s,fs) matrix
    i_loc = sub2ind([N.s,N.k],kron(ones(N.k,1),(1:N.s)'),pol.i_kpr(:));
    G_k = sparse(i_loc,1:n_states,pol.pr_stay(:),n_states,n_states);
    
    % I(k',s | s) entry 
    i_kentry = pol.i_kpr_up(:);
    i_entry_loc = sub2ind([N.s,N.k],(1:N.s)',i_kentry(:));
    G_entry = sparse(i_entry_loc,1:N.s,pol.pr_fstar_enter(:),N.s*N.k,N.s);
%     G_entry = sparse(i_entry_loc,1:N.s,pol.i_entry(:),N.s*N.k,N.s);

    % H(x'|x) matrix (do transpose here)
    H_s = gr.Ps';
    
    % invariant distribution of s
    ps_dist = gr.ps_dist;
    
    % invariant dist of entrants in next period
    g_enter = gr.Ps'*reshape(G_entry*ps_dist(:),N.s,N.k);

    % initialize and get s.s. distribution
    diff_g = 1;
    tol_g = 10^-9;
    maxiter_g = 1000;
    iter_g = 0;
    g = ones(n_states,1)/n_states;
    while diff_g > tol_g && iter_g<maxiter_g

        iter_g = iter_g + 1;

        % move continuers forward
        g_cont = H_s*reshape(G_k*g,N.s,N.k);
        
        % add in entrants
        g1 = g_enter(:) + g_cont(:);

        diff_g = max(abs(g1(:)-g(:)));

        g = g1(:);
    end

    g_ss = reshape(g,N.s,N.k);

end

%% VFI solver (one step induction)
function [Evf,pol] = solve_vfi_onestep(par,N,grids,initvals)

    % unpack parameters
    v2struct(par)
    v2struct(grids)

    % set effective beta
    if par.fixedR==0
        % endogenous R ("GE")
        betaa_pr = betaa;
    else
        % fixed R ("PE")
        betaa_pr = betaa*(P/P_pr);
    end
    
    % next period value function
    Ev = initvals.Ev;
    
    % scrappage recovery function
    kapp  = par.kapp + (1-par.kapp)*(zet==0 && lam==0);
    hk = @(ks) (1-zet)*(1-lam)*(1+ks).^(kapp-1);
    
    %%%%%%% Solve backward induction
    % no resale frictions: directly use analytical solution to compute v
    v_cont = -Q*grids.kstar_noirr_pr + betaa_pr*Ev(i_s_k);

    % pull out value of staying, gross of fixed operating cost
    % (v_stay(s,k))
    v_stay = y + Q*(1-delta).*k_xgrid + repmat(v_cont(:),1,N.k);
    v_up_entry = v_cont; % continuation value for entering
        
    % value of exiting
    v_exit = y + hk(k_xgrid)*Q*(1-delta).*k_xgrid;
    % implied f that makes v_stay>v_exit
    f_star = max(sqrt(eps),v_stay - v_exit);
        % Pr(f'<fstar(s',k'))
    [pr_fstar,Ef_pr] =  cdfu(f_star,mu_fstay_a,mu_fstay_b,npow,eta,s_sh,s_xgrid);
        % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
    lhs = (v_stay - Ef_pr).*pr_fstar;
        % (pi(s',k') + Q*(1-lam)*(1-delta)*k')*(1-Pr(f'<fstar(s',k')))
    rhs = v_exit.*(1-pr_fstar);
        % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
    Ev = Ps*(lhs+rhs);

    % value of entry
    f_star_enter = v_up_entry;
    [pr_fstar_enter,Ef_enter] =  cdfu(f_star_enter,mu_fenter_a,mu_fenter_b,npow,eta,s_sh,s_grid.^((epsi-1)/epsi));

    %%% pull out all policy functions            
    pol.i_kpr = repmat(i_kpr(:),1,N.k);
    pol.i_kpr_up = i_kpr;
    pol.kpr = repmat(grids.kstar_noirr_pr(:),1,N.k);
    pol.kpr_up = grids.kstar_noirr_pr;
    % Policy functions common to both lam=0 and lam>0
    pol.grk = log(pol.kpr(:)./k_xgrid(:));
    pol.fstar = f_star;
    pol.f_star_enter = f_star_enter;
    pol.pr_stay = pr_fstar; % transition probs for stayers
    pol.pr_fstar_enter = pr_fstar_enter; % entry prop conditioned on z draw
    pol.Ef_pr = Ef_pr;
    pol.Ef_enter = Ef_enter;
        % value functions
    pol.v_exit = v_exit;
    pol.v_stay = v_stay;
    Evf = Ev;
    
end

%% Pull out transition matrix / move distribution forward
function [g1,G_k,H_s] = get_trmat(N,gr,pol,trmat,g0)

    %%% check if transition matrices already supplied
    try
        G_k = trmat.G_k;
        H_s = trmat.H_s;
        
    catch
        
        % total # of states
        n_states = N.k*N.s;

        % I(k',s,fs | k,s,fs) matrix
        i_loc = sub2ind([N.s,N.k],kron(ones(N.k,1),(1:N.s)'),pol.i_kpr(:));
        G_k = sparse(i_loc,1:n_states,pol.pr_stay(:),n_states,n_states);

        % I(k',s | s) entry 
        i_kentry = pol.i_kpr_up(:);
        i_entry_loc = sub2ind([N.s,N.k],(1:N.s)',i_kentry(:));
        G_entry = sparse(i_entry_loc,1:N.s,pol.pr_fstar_enter(:),N.s*N.k,N.s);
        
        % invariant distribution of s
        ps_dist = gr.ps_dist;

        % invariant dist of entrants in next period
        g_enter = gr.Ps'*reshape(G_entry*ps_dist(:),N.s,N.k);
        
        % H(x'|x) matrix (do transpose here)
        H_s = gr.Ps';

    end

    %%% iterate forward
    % move continuers forward
    g_cont = H_s*reshape(G_k*g0,N.s,N.k);
    % add in entrants
    g1 = g_enter(:) + g_cont(:);
    

end

%%  Pr(f<fbar) and E[f|f<fbar] functions

% pareto specification (nesting uniform)
function [prx , Ex] = cdfu(f,a,fm,k,sk,sx_sh,sx)
% a = lower-bound, set to 0 by assumption
% fm = upper-bound of pareto
% k = shape parameter of pareto (=1 --> uniform)
% sk = shape parameter cor(f,s)
% sx_sh = location of cor(f,s), (=1 --> default)
% grid of s

    % G(x,s)
    xm = fm.*(sx_sh - (1+sx).^sk(1) + (1+sx).^(sk(1)+sk(2)));
    % CDF 
    f_Prx = @(x) (x./xm).^k;
    % E[x|x<xbar]*Pr(x<xbar)
    f_Ex = @(x) (k/(k+1))./(xm.^k).*(x.^(k+1));

     % force lower bound to be zero
    lb_a = max(a,0);
    
    % compute prob
    prx = f_Prx(f);
    prx(f<lb_a) = 0;
    prx(f>xm) = 1;
    
    % compute conditional expected value
    Ex = f_Ex(f)./prx;
    Ex(f<lb_a) = 0;
    temp = f_Ex(xm);
    Ex(f>xm) = temp(f>xm);
    
end