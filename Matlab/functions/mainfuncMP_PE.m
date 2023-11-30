function output = mainfuncMP_PE(par,N,grids,pol,sw,initval)

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

end

%% Set grids
function grids = setgrid(par,N)

    v2struct(par)
    
    %%% demand externality ("today")
    if chi>0
        Yd = 1/(P*chi);
        Afac = A*P*Yd^(1/epsi);
    else
        % PC and Yd compiled into par structure
        Afac = A*PC*Yd^(1/epsi-1);
    end

    %%% demand externality ("tomorrow")
        % for MIT shocks
        % Note: Afac_pr = Afac if steady-state!
    try
        if chi>0
            Yd_pr = 1/(P_pr*chi);
            Afac_pr = A*P_pr*Yd_pr^(1/epsi);
        else
            disp('Fixed L model not available')
        end
    catch
        Yd_pr = Yd;
        Afac_pr = Afac;
    end
    
    %%%% Get corrected beta's and Q(t-1)
    try
        % Doing MIT shock
        [betaa_pr,betaa_last,Q_last] = get_betaa(betaa,P_last,P,P_pr,Q,chi,par.fixedR);
    catch
        % Doing steady-state
        betaa_pr = betaa;
        betaa_last = betaa;
        Q_last = Q;
    end

    %%% s shocks
    [logs_grid,grids.Ps] = tauchen(N.s,mu_shift,rho_s,sigma_us,2.5);
%     [logs_grid,grids.Ps] = tauchen(N.s,mu_shift,rho_s,sigma_us,3);
    % [Ps, logs_grid] = rouwen(rho_s,mu_shift ,sigma_s , N.s);
    grids.s_grid = exp(logs_grid)*A;
    
%     %%% i.i.d. fixed cost shocks
%     [logfstay_grid,grids.Pfstay] = tauchen(N.fstay,mu_fstay,rho_fstay,sigma_ufstay,3);
%     grids.fstay_grid = exp(logfstay_grid);
%     grids.Pfstay = grids.Pfstay(1,:);
%     
%     %%% s x fixed cost
%     grids.Ps_fstay = kron(grids.Pfstay,grids.Ps);
    
    %%% invariant distribution of s
    ps_dist = grids.Ps^1000;
    grids.ps_dist = par.Me*ps_dist(1,:);
    
    %%% k grids
    if lam>0 || zet>0
        % general case, need to solve for kpr using bellmans
        r = 1/betaa-1;
        try
            k_max = par.kstar;
        catch
    %         kstar = 1050 + (grids.s_grid(end)*alpha_k/(r+delta))^(1/(1-alpha_k));
            k_max = 600;
        end
        
        % k grid
        del_grid = 10;
        grids.k_grid = gen_k_grid(k_max, delta, del_grid, N.k);
        
        % inaction grid
        grids.i_inact = max(1,(1:1:N.k)'- del_grid);

        %%% full mesh grid (needed to compute net y_cont)
        % primary grids
        [grids.s_xgrid0, grids.k_xgrid] =ndgrid(grids.s_grid,grids.k_grid);
        grids.s_xgrid = grids.s_xgrid0.^((epsi-1)/epsi);
        % inaction region
        [~, i_inact] =ndgrid(grids.s_grid,grids.i_inact);
        i_inact = sub2ind([N.s,N.k],repmat((1:N.s)',N.k,1),i_inact(:));
        grids.i_inact_xgrid = reshape(i_inact,N.s,N.k);
        grids.k_inact_xgrid = max((1-delta)*grids.k_xgrid,grids.k_grid(1));
        
    elseif lam==0 && zet==0
        % special case, no need to iterate on bellman
        % Note: for ss, k(t+1) and k(t) are the same grid
        
        %%% set discounting for analytical solution

        
        % override original grid size, replace directly with "on-grid" size
        N.k = N.s;
        
        % compute E[s'|s]
        ss = grids.s_grid.^(((epsi-1)/epsi)*(1+alpha_l/(1-alpha_l)));
        Es = grids.Ps*ss;
        
        % directly obtain the optimal kpr under no irreversibility (t+1)
        vartheta_pr = Q*(1/betaa_pr-1+delta)/alpha_k * Afac_pr.^-(1+alpha_l/(1-alpha_l)) * alpha_l^-(alpha_l/(1-alpha_l));
        grids.kstar_noirr_pr = (vartheta_pr./Es).^(1/(alpha_k-1+alpha_k*alpha_l/(1-alpha_l)));
        % directly back out grid number (trivially 1:Nk)
        grids.i_kpr = 1:N.k;
        grids.i_s_k = sub2ind([N.s,N.k],(1:N.s)',grids.i_kpr(:));
    
        % kpr(t): Use this to construct the "today" grid
        vartheta = Q_last*(1/betaa_last-1+delta)/alpha_k * Afac.^-(1+alpha_l/(1-alpha_l)) * alpha_l^-(alpha_l/(1-alpha_l));
        grids.kstar_noirr = (vartheta./Es).^(1/(alpha_k-1+alpha_k*alpha_l/(1-alpha_l)));
        % override original grid, and replace directly with "on-grid" soln
        grids.k_grid = grids.kstar_noirr;
        %%% for transitional dynamics, override with ss grid for k at t=0
        try
            if t==1
                grids.k_grid = par.k_grid;
                disp('t=1')
            end
        catch
%             disp('solve stationary model')
        end

        %%% full mesh grid (needed to compute net y_cont)
        % primary grids
        [grids.s_xgrid0, grids.k_xgrid] =ndgrid(grids.s_grid,grids.k_grid);
        grids.s_xgrid = grids.s_xgrid0.^((epsi-1)/epsi);

    end
    
    %%% optimal static capacity utilization
    if c==0
        grids.ustar = 1;
    else
        % unconstrained
        xx = (alpha_k.*grids.s_xgrid.*grids.k_xgrid.^(alpha_k-1)/c).^(1/(1-alpha_k)).* (alpha_l.*grids.s_xgrid.*grids.k_xgrid.^alpha_k./par.w).^((1/(1-alpha_l))*(alpha_l/(1-alpha_k)));
        u_uncons = xx.^( (1-alpha_l)*(1-alpha_k) / ((1-alpha_l)*(1-alpha_k)-alpha_l*alpha_k) );
        grids.ustar = min(u_uncons,1);   % cannot exceed 1 
    end
    % optimal labor, given utilization constraint
    grids.lstar = (alpha_l*Afac.*grids.s_xgrid.*grids.k_xgrid.^alpha_k).^(1/(1-alpha_l));    
    % optimal profits (before investment)
    grids.y = Afac.*grids.s_xgrid.*grids.k_xgrid.^alpha_k.*grids.lstar.^alpha_l ...
                - c*grids.ustar.*grids.k_xgrid - grids.lstar;
    % nominal revenue
    grids.py_i = Afac.*grids.s_xgrid.*grids.k_xgrid.^alpha_k.*grids.lstar.^alpha_l;
    % nominal output yi
    grids.y_i = grids.s_xgrid0.*grids.k_xgrid.^alph.*grids.lstar.^(1-alph);
    
    % actual utilized capital
    grids.x_xgrid = grids.ustar.*grids.k_xgrid;

    % reshape into exo-endo
    grids.y = reshape(grids.y,N.s,N.k);
%     grids.y_cont = reshape(grids.y_cont,N.s,N.k);
    grids.py_i = reshape(grids.py_i,N.s,N.k);
    grids.y_i = reshape(grids.y_i,N.s,N.k);
    grids.k_xgrid = reshape(grids.k_xgrid,N.s,N.k);
    grids.x_xgrid = reshape(grids.x_xgrid,N.s,N.k);
    if lam>0 || zet>0
        grids.k_inact_xgrid = reshape(grids.k_inact_xgrid,N.s,N.k);
    end
    
end

%%% Set beta's (function is used only in setgrid.m)
function [betaa_pr,betaa_last,Q_last] = get_betaa(betaa,P_last,P,P_pr,Q,chi,sw_fixedR,PQ_sw)
    % needed for analytical solution for q=Q special case

    if sw_fixedR == 0
        % endogenous R ("GE")

        betaa_last = betaa;
        betaa_pr = betaa;

    elseif sw_fixedR == 1
        % fixed R ("PE")

        % Discounting ("today")
        if chi>0
            % endogenous labor
%             betaa_pr = betaa*(P_pr/P);
%             betaa_pr = betaa*(P/P_pr);
            betaa_pr = betaa;
        else
            % fixed labor
            disp('code does not allow for fixed labor supply')
        end

        % Discounting ("yesterday")
        if chi>0
            % endogenous labor
%             betaa_last = betaa*(P/P_last);
%             betaa_last = betaa*(P_last/P);
            betaa_last = betaa;
        else
            % fixed labor
            disp('code does not allow for fixed labor supply')
        end

    end

    % for MIT shocks ("yesterday's capital price")
    if nargin==8
        if PQ_sw == 1
            % model where P=Q across the transition path
            Q_last = P_last;
        else
            Q_last = Q;
        end
    else
        Q_last = Q;
    end

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

        if lam~=0
            
            % solve up
            i_kpr_up = zeros(N.s,1);
            v_up_cont = zeros(N.s,1);
            % reset: kprime
            last_kpr = 1;
            for is=1:N.s
                % reset: v
                vlast = -1e22;
                for ikpr = last_kpr:N.k

                    % compute candidate value
                    v_candidate = -Q*k_grid(ikpr)/P0 + betaa*Ev(is,ikpr);

                    % compare
                    if v_candidate>=vlast && ikpr<N.k
                        vlast = v_candidate;
                    else
                        v_up_cont(is) = v_candidate;
                        i_kpr_up(is) = ikpr;
                        last_kpr = max(ikpr-1,1);
                        break
                    end

                end
            end
            kpr_up = k_grid(i_kpr_up);
            
            % solve down
            i_kpr_down = zeros(N.s,1);
            v_down_cont = zeros(N.s,1);
            % reset: kprime
            last_kpr = 1;
            for is=1:N.s
                % reset: v
                vlast = -1e22;
                for ikpr = last_kpr:N.k

                    % compute candidate value
                    v_candidate = -Q*(1-lam)*k_grid(ikpr)/P0 + betaa*Ev(is,ikpr);

                    % compare
                    if v_candidate>=vlast && ikpr<N.k
                        vlast = v_candidate;
                    else
                        v_down_cont(is) = v_candidate;
                        i_kpr_down(is) = ikpr;
                        last_kpr = max(ikpr-1,1);
                        break
                    end

                end
            end
            kpr_down = k_grid(i_kpr_down);
            
            
            [~, i_inact_s] =ndgrid(grids.s_grid,grids.i_inact);
            i_inact_s = reshape(sub2ind([N.s,N.k],repmat((1:N.s)',N.k,1),i_inact_s(:)),N.s,N.k);

             %%% check consistency
            % up
            i_up = bsxfun(@gt,kpr_up(:), k_inact_xgrid); 
            v_up_entry = v_up_cont; % continuation value for entering
            v_up_cont = bsxfun(@times,i_up,v_up_cont(:)) + ~i_up*(0);
            % down
            i_down = bsxfun(@lt,kpr_down(:),k_inact_xgrid);
            v_down_cont = bsxfun(@times,i_down,v_down_cont(:)) + ~i_down*(0);
            % inaction
            I_inact = ~i_up & ~i_down;
            v_inact_cont = squeeze(I_inact).*betaa.*Ev(i_inact_s);

            % pull out actual v_cont
            v_cont = v_up_cont + v_down_cont + v_inact_cont;
            
            % pull out value of staying, gross of fixed operating cost
            % (v_stay(s,k))
            v_stay = (y + i_up.*Q*(1-delta).*k_xgrid + i_down.*Q*(1-lam)*(1-delta).*k_xgrid)/P0 + v_cont;

        elseif lam==0 && zet>0
            % no downsizing frictions, but exit frictions exists
            i_kpr = zeros(N.s,1);
            v_cont = zeros(N.s,1);
            last_kpr = 1;
            for is=1:N.s
                % reset: v
                vlast = -1e22;
                for ikpr = last_kpr:N.k

                    % compute candidate value
                    v_candidate = -Q*k_grid(ikpr) + betaa*Ev(is,ikpr);

                    % compare
                    if v_candidate>=vlast && ikpr<N.k
                        vlast = v_candidate;
                    else
                        v_cont(is) = v_candidate;
                        i_kpr(is) = ikpr;
                        last_kpr = max(ikpr-1,1);
                        break
                    end

                end
            end
            
            % pull out value of staying, gross of fixed operating cost
            % (v_stay(s,k))
            v_stay = y + Q*(1-delta).*k_xgrid + repmat(v_cont(:),1,N.k);
            v_up_entry = v_cont; % continuation value for entering
            
        else
            % no resale frictions: directly use analytical solution to compute v
            v_cont = -Q*grids.kstar_noirr + betaa*Ev(i_s_k);
            
            % pull out value of staying, gross of fixed operating cost
            % (v_stay(s,k))
            v_stay = y + Q*(1-delta).*k_xgrid + repmat(v_cont(:),1,N.k);
            v_up_entry = v_cont; % continuation value for entering
            
        end
        
        if ~isnan(mu_fstay_b)
            % model with entry / exit
            
            %%% value of exiting
            v_exit = (y + hk(k_xgrid)*Q*(1-delta).*k_xgrid)/P0;

            %%% implied f that makes v_stay>v_exit
            f_star = max(sqrt(eps),v_stay - v_exit)*P0;
    % f_star=1000000000; % debug
            % Pr(f'<fstar(s',k'))
            [pr_fstar,Ef_pr] =  cdfu(f_star,mu_fstay_a,mu_fstay_b,npow,eta,s_sh,s_xgrid);
            % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
            lhs = (v_stay - Ef_pr/P0).*pr_fstar;
            % (pi(s',k') + Q*(1-lam)*(1-delta)*k')*(1-Pr(f'<fstar(s',k')))
            rhs = v_exit.*(1-pr_fstar);
            % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
            Ev = Ps*(lhs+rhs);
            
            if min(pr_fstar(:)<0)
                disp('error')
                break
            end
            
        else
            
            %%% value of exiting
            v_exit = zeros(size(y));

            %%% implied f that makes v_stay>v_exit
            f_star = zeros(size(v_exit));
    
            % Pr(f'<fstar(s',k'))
            pr_fstar = ones(size(f_star));

            %%% update Ev (nb: v(s,k) -> v(s',k') -> E[v(s',k')|s])
            % E[f'|f'<fstar(s',k')]
            Ef_pr = zeros(size(f_star));
            % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
            lhs = (v_stay - Ef_pr).*pr_fstar;
            % (pi(s',k') + Q*(1-lam)*(1-delta)*k')*(1-Pr(f'<fstar(s',k')))
            rhs = v_exit.*(1-pr_fstar);
            % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
            Ev = Ps*(lhs+rhs);
            
        end
        
        

        % compute difference
        diff_v = max(abs(Ev(:)-Ev0(:)));
        
    end
    
    if ~isnan(mu_fstay_b)
        
        %%% value of entry
        f_star_enter = v_up_entry*P0;
        [pr_fstar_enter,Ef_enter] =  cdfu(f_star_enter,mu_fenter_a,mu_fenter_b,npow,eta,s_sh,s_grid.^((epsi-1)/epsi));
        
    else
        
        %%% value of entry
        f_star_enter = v_up_entry;
        pr_fstar_enter = zeros(size(f_star_enter));
        Ef_enter = zeros(size(f_star_enter)); % expected entry costs
        
    end
    
    %%% pull out all policy functions
    if lam>0 && zet>0
        i_kpr = bsxfun(@times,i_kpr_up(:),i_up) + bsxfun(@times,i_kpr_down(:),i_down) + repmat(i_inact(:)',N.s,1).*I_inact; % pull out actual kpr
        pol.i_kpr = i_kpr;
        pol.i_kpr_up = i_kpr_up;
        pol.kpr = k_grid(i_kpr);
        pol.kpr_up = kpr_up;
        pol.kpr_down = kpr_down;

%         % pull out total dividends
%         pol.d_stay = y_cont + i_up.*Q.*(1-delta).*k_xgrid + i_down.*Q.*(1-lam).*(1-delta).*k_xgrid - i_up.*P.*pol.kpr - Q.*(1-lam).*i_down.*pol.kpr;
%         pol.d_exit = v_exit;
        
    elseif lam==0 && zet>0
        pol.i_kpr = repmat(i_kpr(:),1,N.k);
        pol.i_kpr_up = i_kpr;
        pol.kpr = k_grid(pol.i_kpr);
        pol.kpr_up = k_grid(pol.i_kpr_up);
        
    elseif lam==0 && zet==0
        pol.i_kpr = repmat(i_kpr(:),1,N.k);
        pol.i_kpr_up = i_kpr;
        pol.kpr = repmat(grids.kstar_noirr(:),1,N.k);
        pol.kpr_up = grids.kstar_noirr;

    end
    % Policy functions common to both lam=0 and lam>0
    pol.grk = log(pol.kpr(:)./k_xgrid(:));
    pol.fstar = f_star;
    pol.f_star_enter = f_star_enter;
    pol.pr_stay = pr_fstar; % transition probs for stayers
    pol.pr_fstar_enter = pr_fstar_enter; % entry prop conditioned on z draw
    pol.Ef_pr = Ef_pr;
    pol.Ef_enter = Ef_enter;

    pol.v_exit = v_exit;
    pol.v_stay = v_stay;
    
    Evf = Ev;    
    
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
    
    % next period value function
    Ev = initvals.Ev;
    
    % scrappage recovery function
    kapp  = par.kapp + (1-par.kapp)*(zet==0 && lam==0);
    hk = @(ks) (1-zet)*(1-lam)*(1+ks).^(kapp-1);
    
    %%%%%%% Solve backward induction

    if lam~=0

        % solve up
        i_kpr_up = zeros(N.s,1);
        v_up_cont = zeros(N.s,1);
        % reset: kprime
        last_kpr = 1;
        for is=1:N.s
            % reset: v
            vlast = -1e22;
            for ikpr = last_kpr:N.k

                % compute candidate value
                v_candidate = -Q*k_grid(ikpr)/P + betaa*Ev(is,ikpr);

                % compare
                if v_candidate>=vlast && ikpr<N.k
                    vlast = v_candidate;
                else
                    v_up_cont(is) = v_candidate;
                    i_kpr_up(is) = ikpr;
                    last_kpr = max(ikpr-1,1);
                    break
                end

            end
        end
        kpr_up = k_grid(i_kpr_up);
            
        % solve down
        i_kpr_down = zeros(N.s,1);
        v_down_cont = zeros(N.s,1);
        % reset: kprime
        last_kpr = 1;
        for is=1:N.s
            % reset: v
            vlast = -1e22;
            for ikpr = last_kpr:N.k

                % compute candidate value
                v_candidate = -Q*(1-lam)*k_grid(ikpr)/P + betaa*Ev(is,ikpr);

                % compare
                if v_candidate>=vlast && ikpr<N.k
                    vlast = v_candidate;
                else
                    v_down_cont(is) = v_candidate;
                    i_kpr_down(is) = ikpr;
                    last_kpr = max(ikpr-1,1);
                    break
                end

            end
        end
        kpr_down = k_grid(i_kpr_down);

        [~, i_inact_s] =ndgrid(grids.s_grid,grids.i_inact);
        i_inact_s = reshape(sub2ind([N.s,N.k],repmat((1:N.s)',N.k,1),i_inact_s(:)),N.s,N.k);

         %%% check consistency
        % up
        i_up = bsxfun(@gt,kpr_up(:), k_inact_xgrid); 
        v_up_entry = v_up_cont; % continuation value for entering
        v_up_cont = bsxfun(@times,i_up,v_up_cont(:)) + ~i_up*(0);
        % down
        i_down = bsxfun(@lt,kpr_down(:),k_inact_xgrid);
        v_down_cont = bsxfun(@times,i_down,v_down_cont(:)) + ~i_down*(0);
        % inaction
        I_inact = ~i_up & ~i_down;
        v_inact_cont = squeeze(I_inact).*betaa.*Ev(i_inact_s);

        % pull out actual v_cont
        v_cont = v_up_cont + v_down_cont + v_inact_cont;

        % pull out value of staying, gross of fixed operating cost
        % (v_stay(s,k))
        v_stay = (y + i_up.*Q*(1-delta).*k_xgrid + i_down.*Q*(1-lam)*(1-delta).*k_xgrid)/P + v_cont;
            
    elseif lam==0 && zet>0     
        % no resale frictions
        i_kpr = zeros(N.s,1);
        v_cont = zeros(N.s,1);
        last_kpr = 1;
        for is=1:N.s
            % reset: v
            vlast = -1e22;
            for ikpr = last_kpr:N.k

                % compute candidate value
                v_candidate = -Q*k_grid(ikpr) + betaa*Ev(is,ikpr);

                % compare
                if v_candidate>=vlast && ikpr<N.k
                    vlast = v_candidate;
                else
                    v_cont(is) = v_candidate;
                    i_kpr(is) = ikpr;
                    last_kpr = max(ikpr-1,1);
                    break
                end

            end
        end
        v_stay = y + Q*(1-delta).*k_xgrid + repmat(v_cont(:),1,N.k);
        v_up_entry = v_cont; % continuation value for entering
        
    else
        % no resale frictions: directly use analytical solution to compute v
        v_cont = -Q*grids.kstar_noirr_pr + betaa*Ev(i_s_k);

        % pull out value of staying, gross of fixed operating cost
        % (v_stay(s,k))
        v_stay = y + Q*(1-delta).*k_xgrid + repmat(v_cont(:),1,N.k);
        v_up_entry = v_cont; % continuation value for entering
        
    end
        
    if ~isnan(mu_fstay_b)
        %%% value of exiting
        v_exit = (y + hk(k_xgrid)*Q*(1-delta).*k_xgrid)/P;

        %%% implied f that makes v_stay>v_exit
        f_star = max(sqrt(eps),v_stay - v_exit)*P0;

        % Pr(f'<fstar(s',k'))
        [pr_fstar,Ef_pr] =  cdfu(f_star,mu_fstay_a,mu_fstay_b,npow,eta,s_sh,s_xgrid);
        % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
        lhs = (v_stay - Ef_pr*(1/P0)).*pr_fstar;
        % (pi(s',k') + Q*(1-lam)*(1-delta)*k')*(1-Pr(f'<fstar(s',k')))
        rhs = v_exit.*(1-pr_fstar);
        % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
        Ev = Ps*(lhs+rhs);

        %%% value of entry
        f_star_enter = v_up_entry*P0;
        [pr_fstar_enter,Ef_enter] =  cdfu(f_star_enter,mu_fenter_a,mu_fenter_b,npow,eta,s_sh,s_grid.^((epsi-1)/epsi));
        
    else
        %%% value of exiting
        v_exit = zeros(size(y));

        %%% implied f that makes v_stay>v_exit
        f_star = zeros(size(y));

        % Pr(f'<fstar(s',k'))
        pr_fstar = ones(size(y));

        %%% Note: Here, v_stay and v_exit are computed in for time t values. So
        %%% the step here computes Ev(t), which is extracted and used for
        %%% solving the t-1 problem.
        % E[f'|f'<fstar(s',k')]
        Ef_pr = zeros(size(y));
        % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
        lhs = (v_stay - Ef_pr).*pr_fstar;
        % (pi(s',k') + Q*(1-lam)*(1-delta)*k')*(1-Pr(f'<fstar(s',k')))
        rhs = v_exit.*(1-pr_fstar);
        % (v(s',k') - E[f'|f<fstar(s',k'))*Pr(f'<fstar(s',k'))
        Ev = Ps*(lhs+rhs);

        %%% value of entry
        f_star_enter = v_up_entry;
        pr_fstar_enter = zeros(size(v_up_entry));
        Ef_enter = zeros(size(v_up_entry)); % expected entry costs
        
    end

    %%% pull out all policy functions            
    if lam ~= 0
        
        i_kpr = bsxfun(@times,i_kpr_up(:),i_up) + bsxfun(@times,i_kpr_down(:),i_down) + repmat(i_inact(:)',N.s,1).*I_inact; % pull out actual kpr
        pol.i_kpr = i_kpr;
        pol.i_kpr_up = i_kpr_up;
        pol.kpr = k_grid(i_kpr);
        pol.kpr_up = kpr_up;
        pol.kpr_down = kpr_down;

%         % pull out total dividends
%         pol.d_stay = y_cont + i_up.*Q.*(1-delta).*k_xgrid + i_down.*Q.*(1-lam).*(1-delta).*k_xgrid - i_up.*P.*pol.kpr - Q.*(1-lam).*i_down.*pol.kpr;
%         pol.d_exit = v_exit;
        
    else
        pol.i_kpr = repmat(i_kpr(:),1,N.k);
        pol.i_kpr_up = i_kpr;
        pol.kpr = repmat(grids.kstar_noirr_pr(:),1,N.k);
        pol.kpr_up = grids.kstar_noirr_pr;
    end
    
    % Policy functions common to both lam=0 and lam>0
    pol.grk = log(pol.kpr(:)./k_xgrid(:));
    pol.fstar = f_star;
    pol.f_star_enter = f_star_enter;
    pol.pr_stay = pr_fstar; % transition probs for stayers
    pol.pr_fstar_enter = pr_fstar_enter; % entry prop conditioned on z draw
    pol.Ef_pr = Ef_pr;
    pol.Ef_enter = Ef_enter;

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

% % power distribution specification
% function [prx , Ex] = cdfu(f,a,b,n)
%     % f ~ f^n * ((n+1)/b^(n+1)), f in (a,b)
%     % E[f|f<fbar] = int_{a}{fbar} f*f^n*((n+1)/b^(n+1)) =
%     % (n+1)/((n+2)*b^(n+1))*f^(n+2)
%     
%     % normalizing constant
%     norm_m = (b^(n+1)-a^(n+1))/(n+1);
%     % CDF 
%     f_Prx = @(y) (1/(n+1))*y.^(n+1)/norm_m - (1/(n+1))*a.^(n+1)/norm_m;
%     % E[x|x<xbar]*Pr(x<xbar)
%     if a>=0
%         % if a>0, standard uniform dist
%         f_Ex = @(y) (1/(n+2))*y.^(n+2)/norm_m - (1/(n+2))*a.^(n+2)/norm_m;
%     else
%         % if a<0, mass point at 0 for a<=x<0
%         f_Ex = @(y) (1-f_Prx(0))*(1/(n+2))*y.^(n+2)/norm_m;
%     end
%     % force lower bound to be zero
%     lb_a = max(a,0);
%     
%     % compute prob
%     prx = f_Prx(f);
%     prx(f<lb_a) = 0;
%     prx(f>b) = 1;
%     
%     % compute conditional expected value
%     Ex = f_Ex(f)./prx;
%     Ex(f<lb_a) = 0;
%     Ex(f>b) = f_Ex(b);
% 
% end

% % % % normal specification (nesting uniform)
% % % function [prx , Ex] = cdfu(f,a,fm,k,sk,sx_sh,sx)
% % % % a = lower-bound, set to 0 by assumption
% % % % fm = upper-bound of pareto
% % % % k = shape parameter of pareto (=1 --> uniform)
% % % % sk = shape parameter cor(f,s)
% % % % sx_sh = location of cor(f,s), (=1 --> default)
% % % % grid of s
% % % 
% % % sig=1;
% % % f = max(sqrt(eps),f);
% % % 
% % % 
% % %     % G(x,s)
% % %     xm = fm.*(sx_sh - (1+sx).^sk(1) + (1+sx).^(sk(1)+sk(2)));
% % %     % CDF 
% % %     f_Prx = @(x) normcdf(f,xm);
% % %     % E[x|x<xbar]*Pr(x<xbar)
% % %     f_Ex = @(x) (k/(k+1))./(xm.^k).*(x.^(k+1));
% % % 
% % %      % force lower bound to be zero
% % %     lb_a = max(a,0);
% % %     
% % %     % compute prob
% % %     prx = f_Prx(f);
% % % %     prx(f<lb_a) = 0;
% % % %     prx(f>xm) = 1;
% % %     
% % %     
% % %     % compute conditional expected value
% % %     num = (log(f)-xm-sig^2)/sig;
% % %     denom = (log(f)-xm)/sig;
% % %     Ex = exp(xm+sig^2/2).*(normcdf(num,xm)./normcdf(denom,xm));
% % % %     Ex(f<lb_a) = 0;
% % % %     temp = f_Ex(xm);
% % % %     Ex(f>xm) = temp(f>xm);
% % %     
% % % end

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