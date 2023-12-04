% f_ss_covid.m
%
%   This is the basic function to solve the model after Covid shock
%
%   This function solves the model by varying P to clear C. This is exactly
%   the same function as f_ss_trade_shock.m, except that this is used
%   together with f_PE_ss_autarky.m; the outputs are all therefore also
%   tagged with _PE.
%
%   There are three inputs:
%       1. Input parameters
%       2. Location to export mat files (optional)
%       3. Suffix to tag onto file name (optional)

function output_moment = f_PE_ss_covid(par_in,file_out,file_suffix)

%% parameters and grids: set baseline parameter
disp("[f_PE_ss_covid.m] Setting up parameters & grid...")

%%% retrieve baseline parameters and grid sizes
try
    output  = mainfuncMP([],[],[],[],1,par_in);
catch
    disp('no parameters supplied. using default')
    output  = mainfuncMP([],[],[],[],1);
end
par     = output.par;
N       = output.N;
clear output

%%% load capital grid from autarky solution
try
    par.kstar  = par_in.kmax_autarky;
catch
    disp('no new kmax provided');
    par.kstar = 300;
end

%%% load new measure from autarky solution
par.Me = par_in.Me;

%%% load pre-shock prices
Pnow = par_in.P_autarky;
Pss_notrade = par_in.P_autarky;

%%% set parameters for foreign imports
try
    dM = par_in.dM;
    Pf = par_in.Pf;
catch
    disp('no parameters supplied for foreign firms')
    dM = .2; % measure of foreign firms
    Pf = Pss_notrade; % price of foreign goods
end

%% Solve GE

disp("[f_PE_ss_covid.m] We're going to vary P to clear C, holding fix M at autarky value.")
disp("[f_PE_ss_covid.m] Solve GE via bisection...")

% form par2 structure
par2 = par;

%%% Solve GE
P_guess = zeros(1000,1);
new_C_store = zeros(1000,1);
itermax = 1000;
tolw = 1e-4;
relx = 0.01;
for iter=1:itermax

    % set Q to P for now
    par2.P  = Pnow;
    if par2.lam==0 && par2.zet==0
        par2.Q  = par_in.Q;
    else
        par2.Q  = Pss_notrade;
    end

    %%% set the grid
    if par2.lam==0 && par2.zet==0
        gr = mainfuncMP_q1(par2,N,[],[],2);
    else
        gr = mainfuncMP(par2,N,[],[],2);
    end
    % replace Nk if special case
    if par2.lam==0 && par2.zet==0
        N.k = length(gr.k_grid);
    end

    %%% solve the individual's problem
    if par2.lam==0 && par2.zet==0
        % no irreversibility
        try
            output  = mainfuncMP_q1(par2,N,gr,[],3,initval);
        catch
            output  = mainfuncMP_q1(par2,N,gr,[],3);
        end
    else
        try
            output  = mainfuncMP(par2,N,gr,[],3,initval);
        catch
            output  = mainfuncMP(par2,N,gr,[],3);
        end
    end
    pol     = output.pol;
    Ev       = output.Ev;
    clear output

    %%% get distribution
    output  = mainfuncMP(par2,N,gr,pol,4);
    Px      = output;
    clear output

    %%% Compute consumption for old and new
    % implied total consumption given guess of P
    old_C = 1/(par2.P*par2.chi);
    % inner integral for domestic goods
    inner_Cd = gr.y_i(:).^((par2.epsi-1)/par2.epsi)'*Px.g_ss(:);
    % inner integral for foreign goods
    inner_Cf = (Pf.^(-par2.epsi)*(par2.P^par2.epsi)*old_C)^((par2.epsi-1)/par2.epsi)*dM;
    cf_i = Pf.^(-par2.epsi)*(par2.P^par2.epsi)*old_C; % individual foreign firm
    % actual consumption
    new_C = (inner_Cd + inner_Cf)^(par2.epsi/(par2.epsi-1));
    % update value function
    initval.Ev = Ev;

    % check tolerance
    diff_C = new_C - old_C;
    if abs(diff_C) < tolw
        disp(diff_C);
        break
    end

    % update prices
    new_C_store(iter) = new_C;
%     P_guess(iter) = 1/(((1-relx)*old_C + relx*new_C)*par2.chi);
    if iter==1

        lc = min(new_C,old_C);
        rc = max(new_C,old_C);

    else

        lc_cand = min(new_C,old_C);
        rc_cand = max(new_C,old_C);

        lc = max(lc,lc_cand);
        rc = min(rc,rc_cand);

    end
    P_guess(iter) = 1/(((1-relx)*lc + relx*rc)*par2.chi);
    Pnow = P_guess(iter);

    % In stead of plotting with interupt the screen with popups, I just print things out here
    disp(sprintf("     Iteration %d: C = %d, 𝚫C = %d, P = %d", iter, new_C, Pnow, diff_C))

    % % plot progress
    % subplot(311)
    % drawnow
    % hold on
    % plot(iter,diff_C,'x');
    % subplot(312)
    % plot(iter,Pnow,'x');
    % drawnow
    % hold on
    % subplot(313)
    % plot(iter,new_C,'x');
    % drawnow
    % hold on

    % % print to screen progression
    % disp(diff_C);

    if abs(diff_C)<10
        relx=.1;
    elseif abs(diff_C)<1
        relx = .2;
    elseif abs(diff_C)<0.5
        relx = .5;
    end

end
close all

% compute moments
par2 = par;
par2.P  = Pnow;
if par2.lam==0 && par2.zet==0
    par2.Q  = par_in.Q;
else
    par2.Q  = Pss_notrade;
end
par2.Pf = Pf;
par2.dM = dM;
pol.cf_i = cf_i;
output_moment = get_moments(par2,N,gr,Px,pol,1);


% export full workspace
if nargin==3
    % file path and suffix supplied
    if par2.lam>0
        save([file_out '/ss_covid' file_suffix])
        save([file_out '/Pss_covid' file_suffix],'Pnow')
    else
        save([file_out '/ss_covid_q1' file_suffix])
        save([file_out '/Pss_covid_q1' file_suffix],'Pnow')
    end

elseif nargin==2
    % file path supplied
    if par2.lam>0
        save([file_out '/ss_covid'])
        save([file_out '/Pss_covid'],'Pnow')
    else
        save([file_out '/ss_covid_q1'])
        save([file_out '/Pss_covid_q1'],'Pnow')
    end

else
    % no file path supplied
    mkdir('mat_files')
    if par2.lam>0
        save('mat_files/ss_covid')
        save('mat_files/Pss_covid','Pnow')
    else
        save('mat_files/ss_covid_q1')
        save('mat_files/Pss_covid_q1','Pnow')
    end
end



end
