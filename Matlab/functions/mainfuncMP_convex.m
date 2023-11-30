function output = mainfuncMP_convex(par,N,grids,pol,sw,initval)
% This version: Kitchen sink. All cost.

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
    % directly repackage
    par = Nnew;
    
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
    par.c1 = Nnew.c1;
    par.C_f = Nnew.C_f;
    
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

end

%% Set grids
function grids = setgrid(par,N)

    v2struct(par)
    
    % demand externality ("today")
    Yd = 1/(P*chi);
    Afac = A*P*Yd^(1/epsi);

    %%% s shocks
    [logs_grid,grids.Ps] = tauchen(N.s,mu_shift,rho_s,sigma_us,2.5);
    grids.s_grid = exp(logs_grid)*A;
    
    %%% invariant distribution of s
    ps_dist = grids.Ps^1000;
    grids.ps_dist = par.Me*ps_dist(1,:);
    
    %%% k grids
    k_max = par.kstar;
    grids.k_grid = exp(linspace(log(1+1e-4),log(1+k_max),N.k))'-1;

    %%% full mesh grid (needed to compute net y_cont)
    % primary grids
        % s_xgrid0: physical tfp
        % s_xgrid: revenue tfp (omega)
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

       % solve continuers investment problem
        kpr_adj_up = zeros(N.s,N.k);
        kpr_adj_down = zeros(N.s,N.k);
        v_stay_up = zeros(N.s,N.k);
        v_stay_down = zeros(N.s,N.k);
        kpr_0 = zeros(N.s,N.k);
        v_stay_0 = zeros(N.s,N.k);
        for is=1:N.s
            
            % fit interpolating function
            F_Ev = griddedInterpolant(k_grid,Ev(is,:),'makima');
            
            % up
            [kpr_adj_up(is,:), v_stay_up(is,:)] = goldenx(@vv_up,k_grid(1)*ones(size(k_grid)),k_grid(N.k)*ones(size(k_grid)),y(is,:)',F_Ev,k_grid,Q,lam,c1,delta,betaa,C_f);
            
            % down
            [kpr_adj_down(is,:), v_stay_down(is,:)] = goldenx(@vv_down,k_grid(1)*ones(size(k_grid)),k_grid(N.k)*ones(size(k_grid)),y(is,:)',F_Ev,k_grid,Q,lam,c1,delta,betaa,C_f);
            
            % inaction
            kpr_0(is,:) = max(k_grid(1),(1-delta)*k_grid);
            v_stay_0(is,:) = vv(kpr_0(is,:)',y(is,:)',F_Ev,k_grid,Q,0,c1,delta,betaa,0*C_f);
        end
        % check consistency
            % up
        i_adj_up = kpr_adj_up > max(k_grid(1),(1-delta)*k_grid');
        v_stay_up(~i_adj_up) = -1e10;
            % down
        i_adj_down = kpr_adj_down < max(k_grid(1),(1-delta)*k_grid');
        v_stay_down(~i_adj_down) = -1e10;
            % get adjuster's value
        v_stay_adj = max(v_stay_up,v_stay_down);
        i_stay_adj = v_stay_up>=v_stay_down;
            % adjuster's investment
        kpr_adj = i_stay_adj.*kpr_adj_up + ~i_stay_adj.*kpr_adj_down;
            % max adj vs don't adj
        i_val_adj = (v_stay_adj>v_stay_0);
        v_stay = i_val_adj.*v_stay_adj + ~i_val_adj.*v_stay_0;
            % actual investment function
        kpr = i_val_adj.*kpr_adj + ~i_val_adj.*kpr_0;
        
        % solve exit vs stay decision
            % value of exiting
        v_exit = y + hk(k_xgrid)*Q*(1-delta).*k_xgrid - C_f*Q*(1-delta).*k_xgrid - c2*(1-delta)^2*k_xgrid;
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

        % check: probs>=0
        if min(pr_fstar(:)<0)
            disp('error')
            break
        end

        % check value function tol
        diff_v = max(abs(Ev(:)-Ev0(:)));
        
    end
    
    % entry
    kpr_up = zeros(N.s,1);
    v_up_entry = zeros(N.s,1);
    for is=1:N.s
        F_Ev = griddedInterpolant(k_grid,Ev(is,:),'pchip');
        [kpr_up(is), v_up_entry(is)] = goldenx(@vv,k_grid(1),k_grid(N.k),0,F_Ev,sqrt(eps),Q,0,c1*0,delta,betaa,C_f);
    end
        % value of entry indifference condition
    f_star_enter = v_up_entry;
    [pr_fstar_enter,Ef_enter] =  cdfu(f_star_enter,mu_fenter_a,mu_fenter_b,npow,eta,s_sh,s_grid.^((epsi-1)/epsi));
    
    %%% pull out all policy functions
    pol.kpr = kpr;
    pol.kpr_up = kpr_up;
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

% "combined" value function
function v_candidate = vv(kp,y_is,F_Ev_is,k_grid,Q,lam,c1,delta,betaa,C_f)
    inv = kp - (1-delta)*k_grid;
    xi = (exp(log(kp)-log(k_grid)) - (1-delta)).^2.*k_grid;
    Qq = Q*(1-lam*(inv<0));
    v_candidate = y_is - Qq.*inv - c1*xi - Q*(1-delta)*C_f*k_grid + betaa*F_Ev_is(kp);
end

% value of investment
function v_candidate = vv_up(kp,y_is,F_Ev_is,k_grid,Q,lam,c1,delta,betaa,C_f)
    inv = kp - (1-delta)*k_grid;
    xi = (exp(log(kp)-log(k_grid)) - (1-delta)).^2.*k_grid;
%     xi = (kp./k_grid - (1-delta)).^2.*k_grid;
    Qq = Q;
    v_candidate = y_is - Qq.*inv - c1*xi - Q*(1-delta)*C_f*k_grid + betaa*F_Ev_is(kp);
end

% value of disinvestment
function v_candidate = vv_down(kp,y_is,F_Ev_is,k_grid,Q,lam,c1,delta,betaa,C_f)
    inv = kp - (1-delta)*k_grid;
    xi = (exp(log(kp)-log(k_grid)) - (1-delta)).^2.*k_grid;
%     xi = (kp./k_grid - (1-delta)).^2.*k_grid;
%     Qq = Q*(1-lam);
    v_candidate = y_is + Q*(1-lam)*(1-delta)*k_grid - Q*(1-lam).*kp - c1*xi - Q*(1-delta)*C_f*k_grid + betaa*F_Ev_is(kp);
end

%% Solve for s.s. dist
function [g_ss,G_k,H_s] = get_dist(N,gr,pol)

    % total # of states
    n_states = N.k*N.s;
    
    % log grids
    log_k = log(1+gr.k_grid);
    log_kp = log(1+pol.kpr);
    log_kp_up = log(1+pol.kpr_up);
    
    %%% continuers
    %bin capital
    i_left = discretize(log_kp,log_k);
    i_right = min(i_left+1,N.k);
    d_left = log_kp - log_k(i_left);
    d_right = log_k(i_right) - log_kp;
    d_all = log_k(i_right) - log_k(i_left);
    s_left = (d_all-d_left)./d_all;
    s_right = (d_all-d_right)./d_all;
    % I(k',s,fs | k,s,fs) matrix
    i_loc_l = sub2ind([N.s,N.k],kron(ones(N.k,1),(1:N.s)'),i_left(:));
    i_loc_r = sub2ind([N.s,N.k],kron(ones(N.k,1),(1:N.s)'),i_right(:));
    G_k_l = sparse(i_loc_l,1:n_states,s_left(:).*pol.pr_stay(:),n_states,n_states);
    G_k_r = sparse(i_loc_r,1:n_states,s_right(:).*pol.pr_stay(:),n_states,n_states);
    G_k = G_k_l + G_k_r;
    
    clear i_* d_* s_*
    
    %%% entrants
    %bin capital
    i_left = discretize(log_kp_up,log_k);
    i_right = min(i_left+1,N.k);
    d_left = log_kp_up - log_k(i_left);
    d_right = log_k(i_right) - log_kp_up;
    d_all = log_k(i_right) - log_k(i_left);
    s_left = (d_all-d_left)./d_all;
    s_right = (d_all-d_right)./d_all;
    % I(k',s | s) entry 
    i_entry_loc_l = sub2ind([N.s,N.k],(1:N.s)',i_left(:));
    i_entry_loc_r = sub2ind([N.s,N.k],(1:N.s)',i_right(:));
    G_entry = sparse(i_entry_loc_l,1:N.s,s_left(:).*pol.pr_fstar_enter(:),N.s*N.k,N.s) + sparse(i_entry_loc_r,1:N.s,s_right(:).*pol.pr_fstar_enter(:),N.s*N.k,N.s);

    % H(x'|x) matrix (do transpose here)
    H_s = gr.Ps';
    
    % invariant distribution of s
    ps_dist = gr.ps_dist;
    
    % invariant dist of entrants in next period
    g_enter = gr.Ps'*reshape(G_entry*ps_dist(:),N.s,N.k);

    % initialize and get s.s. distribution
    diff_g = 1;
    tol_g = 10^-11;
    maxiter_g = 10000;
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
    
    % solve continuers
    kpr_adj_up = zeros(N.s,N.k);
    kpr_adj_down = zeros(N.s,N.k);
    v_stay_up = zeros(N.s,N.k);
    v_stay_down = zeros(N.s,N.k);
    kpr_0 = zeros(N.s,N.k);
    v_stay_0 = zeros(N.s,N.k);
    for is=1:N.s
        
        % fit interpolating function
        F_Ev = griddedInterpolant(k_grid,Ev(is,:),'makima');
        
        % up
        [kpr_adj_up(is,:), v_stay_up(is,:)] = goldenx(@vv_up,k_grid(1)*ones(size(k_grid)),k_grid(N.k)*ones(size(k_grid)),y(is,:)',F_Ev,k_grid,Q,lam,c1,delta,betaa_pr,C_f);

        % down
        [kpr_adj_down(is,:), v_stay_down(is,:)] = goldenx(@vv_down,k_grid(1)*ones(size(k_grid)),k_grid(N.k)*ones(size(k_grid)),y(is,:)',F_Ev,k_grid,Q,lam,c1,delta,betaa_pr,C_f);
            
        % inaction
        kpr_0(is,:) = max(k_grid(1),(1-delta)*k_grid);
        v_stay_0(is,:) = vv(kpr_0(is,:)',y(is,:)',F_Ev,k_grid,Q,0,c1,delta,betaa_pr,0*C_f);
    end
    
    % check consistency
        % up
    i_adj_up = kpr_adj_up > max(k_grid(1),(1-delta)*k_grid');
    v_stay_up(~i_adj_up) = -1e10;
        % down
    i_adj_down = kpr_adj_down < max(k_grid(1),(1-delta)*k_grid');
    v_stay_down(~i_adj_down) = -1e10;
        % max for adj
    v_stay_adj = max(v_stay_up,v_stay_down);
    i_stay_adj = v_stay_up>=v_stay_down;
        % adj investment
    kpr_adj = i_stay_adj.*kpr_adj_up + ~i_stay_adj.*kpr_adj_down;
    i_val_adj = (v_stay_adj>v_stay_0);
        % vs inaction
    v_stay = i_val_adj.*v_stay_adj + ~i_val_adj.*v_stay_0;
    kpr = i_val_adj.*kpr_adj + ~i_val_adj.*kpr_0;
    
    % value of exiting
    v_exit = y + hk(k_xgrid)*Q*(1-delta).*k_xgrid - c2*(1-delta)^2*k_xgrid;
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

    if min(pr_fstar(:)<0)
        disp('error')
    end

    % entry
    kpr_up = zeros(N.s,1);
    v_up_entry = zeros(N.s,1);
    for is=1:N.s
        F_Ev = griddedInterpolant(k_grid,Ev(is,:),'pchip');
        [kpr_up(is), v_up_entry(is)] = goldenx(@vv,k_grid(1),k_grid(N.k),0,F_Ev,sqrt(eps),Q,0,c1*0,delta,betaa_pr,C_f);
    end
    
    % value of entry
    f_star_enter = v_up_entry;
    [pr_fstar_enter,Ef_enter] =  cdfu(f_star_enter,mu_fenter_a,mu_fenter_b,npow,eta,s_sh,s_grid.^((epsi-1)/epsi));

    %%% pull out all policy functions            
    pol.kpr = kpr;
    pol.kpr_up = kpr_up;
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
        
        % log grids
        log_k = log(1+gr.k_grid);
        log_kp = log(1+pol.kpr);
        log_kp_up = log(1+pol.kpr_up);

        %%% continuers
        %bin capital
        i_left = discretize(log_kp,log_k);
        i_right = min(i_left+1,N.k);
        d_left = log_kp - log_k(i_left);
        d_right = log_k(i_right) - log_kp;
        d_all = log_k(i_right) - log_k(i_left);
        s_left = (d_all-d_left)./d_all;
        s_right = (d_all-d_right)./d_all;
        % I(k',s,fs | k,s,fs) matrix
        i_loc_l = sub2ind([N.s,N.k],kron(ones(N.k,1),(1:N.s)'),i_left(:));
        i_loc_r = sub2ind([N.s,N.k],kron(ones(N.k,1),(1:N.s)'),i_right(:));
        G_k_l = sparse(i_loc_l,1:n_states,s_left(:).*pol.pr_stay(:),n_states,n_states);
        G_k_r = sparse(i_loc_r,1:n_states,s_right(:).*pol.pr_stay(:),n_states,n_states);
        G_k = G_k_l + G_k_r;

        clear i_* d_* s_*

        %%% entrants
        %bin capital
        i_left = discretize(log_kp_up,log_k);
        i_right = min(i_left+1,N.k);
        d_left = log_kp_up - log_k(i_left);
        d_right = log_k(i_right) - log_kp_up;
        d_all = log_k(i_right) - log_k(i_left);
        s_left = (d_all-d_left)./d_all;
        s_right = (d_all-d_right)./d_all;
        % I(k',s | s) entry 
        i_entry_loc_l = sub2ind([N.s,N.k],(1:N.s)',i_left(:));
        i_entry_loc_r = sub2ind([N.s,N.k],(1:N.s)',i_right(:));
        G_entry = sparse(i_entry_loc_l,1:N.s,s_left(:).*pol.pr_fstar_enter(:),N.s*N.k,N.s) + sparse(i_entry_loc_r,1:N.s,s_right(:).*pol.pr_fstar_enter(:),N.s*N.k,N.s);

        % H(x'|x) matrix (do transpose here)
        H_s = gr.Ps';

        % invariant distribution of s
        ps_dist = gr.ps_dist;

        % invariant dist of entrants in next period
        g_enter = gr.Ps'*reshape(G_entry*ps_dist(:),N.s,N.k);

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