%% thresholds.m
%
%   Use this to compute and plot the thresholds and distributions
%
%   Given any policy function and grid, compute the 
%       1. Exit thresholds (for given fx)
%       2. Inaction regions

function out=thresholds(par, N, pol_new, gr, g_ss, closed_moments, sudden_open_moments, fx, pr0_x)


%% Compute exit thresholds for a given fixed cost

% grid the fstar function
F_fx_q1_0 = griddedInterpolant({gr.s_grid,gr.k_grid},pol_new.fstar);

% construct a loss function
yout_q1_0 =@(si) -1*(fx-F_fx_q1_0(si,gr.k_grid)).^2;

% compute s* given the grid of k
s_fx_star = goldenx(yout_q1_0,gr.s_grid(1)*ones(size(gr.k_grid)),gr.s_grid(end)*ones(size(gr.k_grid)));

%% compute s* given pr(exit)
% compute this using the probit outcomes

try
    pr0 = pr0_x;
catch
    pr0 = 0.90;
end

% closed moments
if exist('closed_moments','var')
    if ~isempty(sudden_open_moments)
        [c,s,k] = ndgrid(1,gr.s_grid,gr.k_grid);
        xr = [c(:) log(s(:)) log(k(:))]*closed_moments.bhat_probit;
        xr = reshape(xr,N.s,N.k);
        pr = normcdf(xr);
        ks = zeros(N.s,1);
        for ii=1:N.s
            fpr = griddedInterpolant(log(gr.k_grid),pr(ii,:),'pchip');
            logktemp_r = log(gr.k_grid(end));
            logktemp_l = log(gr.k_grid(1));    
            pr_min = fpr(log(gr.k_grid(1)));
            pr_max = fpr(log(gr.k_grid(end)));

            if pr_min > pr0
                logktemp = log(gr.k_grid(1));
            elseif pr_max < pr0
                logktemp = log(gr.k_grid(end));
            else
                while abs(logktemp_l-logktemp_r)>1e-6

                    [logktemp,~] = goldenx(@(x) -1*(fpr(x)-pr0)^2,logktemp_l,logktemp_r);
                    if fpr(logktemp) > pr0
                        logktemp_r = logktemp;
                    else
                        logktemp_l = logktemp; 
                    end
                end
            end

            ks(ii) = exp(logktemp);
        end
        out.closed_k_prob = ks;
    end
end

% opened moments
if exist('sudden_open_moments','var')
    if ~isempty(sudden_open_moments)
        [c,s,k] = ndgrid(1,gr.s_grid,gr.k_grid);
        xr = [c(:) log(s(:)) log(k(:))]*sudden_open_moments.bhat_probit;
        xr = reshape(xr,N.s,N.k);
        pr = normcdf(xr);
        ks = zeros(N.s,1);
        for ii=1:N.s
            fpr = griddedInterpolant(log(gr.k_grid),pr(ii,:),'pchip');
            logktemp_r = log(gr.k_grid(end));
            logktemp_l = log(gr.k_grid(1));    
            pr_min = fpr(log(gr.k_grid(1)));
            pr_max = fpr(log(gr.k_grid(end)));

            if pr_min > pr0
                logktemp = log(gr.k_grid(1));
            elseif pr_max < pr0
                logktemp = log(gr.k_grid(end));
            else
                while abs(logktemp_l-logktemp_r)>1e-6

                    [logktemp,~] = goldenx(@(x) -1*(fpr(x)-pr0)^2,logktemp_l,logktemp_r);
                    if fpr(logktemp) > pr0
                        logktemp_r = logktemp;
                    else
                        logktemp_l = logktemp; 
                    end
                end
            end

            ks(ii) = exp(logktemp);
        end
        out.opened_k_prob = ks;
    end
end

%% Compute inaction region

% reshape
kpr_s = reshape(pol_new.kpr,N.s,N.k);
k_inact = max((1-par.delta)*gr.k_grid,gr.k_grid(1));

%%% Solve investment / disinvestment thresholds
k_up_threshold = zeros(N.s,1);
k_down_threshold = zeros(N.s,1);
try
    i_irr = par.lam>0 || par.zet>0;
catch
    i_irr = 0;
end
try
    i_fc = par.C_f>0 || par.C_f_zeta>0;
catch
    i_fc = 0;
end
    
if i_irr || i_fc
    for ii=1:N.s
        diff_k = kpr_s(ii,:)-k_inact(:)';
        diff_k(abs(diff_k)<sqrt(eps)) = 0;
        try
            k_up_threshold(ii) = gr.k_grid(find(diff_k>0,1,'last'));
        catch
            k_up_threshold(ii) = gr.k_grid(find(diff_k>=0,1,'last'));
        end
        try
            i_down = find(diff_k<0,1,'first');
%             i_check = any(diff_k(i_down:end)==0);
%             while i_check || i_down==length(diff_k)
%                 i_check = any(diff_k(i_down:end)==0);
%                 i_down = find(diff_k(i_down+1:end)<0,1,'first');
%             end
            k_down_threshold(ii) = gr.k_grid(i_down);
        catch
            k_down_threshold(ii) = gr.k_grid(end);
        end
    end
    k_up_threshold = sort(k_up_threshold) + linspace(1e-4,2e-4,N.s)';
    k_down_threshold = sort(k_down_threshold) + linspace(1e-4,2e-4,N.s)';

    
    out.k_up_threshold = k_up_threshold;
    out.k_down_threshold = k_down_threshold;
else
    kpr_star = griddedInterpolant(gr.s_grid,kpr_s(:,1),'pchip');
    out.s_k0_threshold = zeros(N.k,1);
    for ik=1:N.k
        diff_k = @(ss) -1*(kpr_star(ss)-k_inact(ik)).^2;
%         out.s_k0_threshold(ik) = max(ik*sqrt(eps),fzero(diff_k,k_inact(ik)));
        out.s_k0_threshold(ik) = max(ik*sqrt(eps),goldenx(diff_k,gr.s_grid(1),gr.s_grid(end)));
        if ik>1
            if out.s_k0_threshold(ik)==out.s_k0_threshold(ik-1)
                out.s_k0_threshold(ik) = out.s_k0_threshold(ik)+sqrt(eps);
            end
        end
    end
end

out.s_fx_star = s_fx_star;

% % % plot(gr.k_grid,out.s_k0_threshold)
% % % f_k_inact = griddedInterpolant(gr.s_grid,k_inact(:,1),'pchip');
% % % s_star = f_k_inact(gr.k_grid);
% % % figure
% % % plot(gr.k_grid,s_star)



%% Regression, and then prediction after regression

try 
    epsi = par.epsi;
    log_tfpr = log(gr.s_grid.^((epsi-1)/epsi));
    [ss,kk] = ndgrid(log_tfpr,log(gr.k_grid));

    pr_fixed = .75;

    pr_close = closed_moments.bet_hat_all'*[ones(size(ss(:)))  ss(:) kk(:)]';
    % pr_close = exp(pr_close)-sqrt(eps);
    pr_close = reshape(pr_close,N.s,N.k);
    Fpr_close = griddedInterpolant({log_tfpr,log(gr.k_grid)},pr_close,'linear');


    Q = @(s_x,pr_in) -1*(Fpr_close(s_Xx,log(gr.k_grid)) - pr_in).^2;
    s_closed = goldenx(Q,ones(size(log(gr.k_grid)))*log_tfpr(1),ones(size(log(gr.k_grid)))*log_tfpr(end),pr_fixed);

    pr_open = sudden_open_moments.bet_hat_all'*[ones(size(ss(:)))  ss(:) kk(:)]';
    % pr_open = exp(pr_open)-sqrt(eps);
    pr_open = reshape(pr_open,N.s,N.k);
    Fpr_open = griddedInterpolant({log_tfpr,log(gr.k_grid)},pr_open,'linear');

    Q = @(s_x,pr_in) -1*(Fpr_open(s_x,log(gr.k_grid)) - pr_in).^2;
    s_open = goldenx(Q,ones(size(log(gr.k_grid)))*log_tfpr(1),ones(size(log(gr.k_grid)))*log_tfpr(end),pr_fixed);

    disp([closed_moments.bet_hat_all(:) sudden_open_moments.bet_hat_all(:)])

    figure
    subplot(211)
    suptitle('Predicted probabilities of survival')
    mesh(kk,ss,pr_close);
    title('Close economy')
    subplot(212)
    mesh(kk,ss,pr_open);
    title('Period 0 opening economy')

    figure
    subplot(211)
    suptitle('Predicted probabilities of survival')
    contour(kk,ss,pr_close);
    title('Close economy')
    subplot(212)
    contour(kk,ss,pr_open);
    title('Period 0 opening economy')

    figure
    hold on
    plot(log(gr.k_grid),(s_closed))
    plot(log(gr.k_grid),(s_open))
    xlabel('log k')
    ylabel('log s');
    contour((kk),(ss),g_sk/sum(g_ss(:)));
    legend({'Closed','Sudden opening','Steady-state distribution'},'location','southwest')
    axis([-20 5 -2 2])
    hold off


catch
    
end




