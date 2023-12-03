% f_PE_ss_trade_shock_convex.m
%
%   This is the basic function to solve the model after opening to trade
%
%   This function solves the model by varying P to clear C, holding fix M
%   at autarky value.
%
%   There are three inputs:
%       1. Input parameters
%       2. Location to export mat files (optional)
%       3. Suffix to tag onto file name (optional)

function output_moment = f_PE_ss_trade_shock_convex(par_in,file_out,file_suffix)

%% parameters and grids: set baseline parameter
disp("[f_PE_ss_trade_shock_convex.m] Setting up parameters & grid...")

%%% retrieve baseline parameters and grid sizes
output  = mainfuncMP_convex([],[],[],[],1,par_in); 
par     = output.par;
N       = output.N;
clear output

%%% load capital grid from autarky solution
par.kstar  = par_in.kmax_autarky;

%%% load new measure from autarky solution
par.Me = par_in.Me;

%%% load pre-shock prices
Pnow = par_in.P_autarky;
Pss_notrade = par_in.P_autarky;

%%% set parameters for foreign imports
    % measure of foreign firms
dM = par_in.dM;
    % price of foreign goods
Pf = par_in.Pf;

%% Solve GE
   
% form par2 structure to avoid changing parameter input structure
par2 = par;

%%% Solve GE via bisection
disp("[f_PE_ss_trade_shock_convex.m] We're going to vary P to clear C, holding fix M at autarky value.")
disp("[f_PE_ss_trade_shock_convex.m] Solve GE via bisection...")

P_guess = zeros(1000,1);
new_C_store = zeros(1000,1);
itermax = 1000;
    % tolerance
tolw = 1e-3;
    % dampening
relx = 0.01;
for iter=1:itermax
    
    % Update P
    par2.P  = Pnow;
    
    % Impose Q = price at autarky
    par2.Q  = Pss_notrade;

    %%% set the grid
    gr = mainfuncMP_convex(par2,N,[],[],2);

    %%% solve the individual's problem
    try
        output  = mainfuncMP_convex(par2,N,gr,[],3,initval);
    catch
        % no initial value function supplied
        output  = mainfuncMP_convex(par2,N,gr,[],3);
    end
    pol     = output.pol;
    Ev       = output.Ev;
    clear output

    %%% get distribution
    output  = mainfuncMP_convex(par2,N,gr,pol,4);
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

    % check market clearing error
    diff_C = new_C - old_C;
    if abs(diff_C) < tolw
        disp(diff_C);
        break
    end
    
    % update prices via bisection
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

    % plot progress
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

    
    % print progress to screen
    % disp(diff_C);
    
    % update tolerance as market clearing error declines
    if abs(diff_C)<10
        relx=.1;
    elseif abs(diff_C)<1
        relx = .2;
    elseif abs(diff_C)<0.5
        relx = .5;
    end
        
end
close all

%% compute moments
disp("[f_PE_ss_trade_shock_convex.m] Computing the relevant moments...")

par2 = par;
par2.P  = Pnow;
par2.Q  = Pss_notrade;
par2.Pf = Pf;
par2.dM = dM;
pol.cf_i = cf_i;
output_moment = get_moments_convex(par2,N,gr,Px,pol,1);


%% export full workspace
disp(sprintf("[f_PE_ss_trade_shock_convex.m] export the workspace to %d/ss_tradeshock or Pss_tradeshock...", file_out))

if nargin==3
    % file path and suffix supplied
    save([file_out '/ss_tradeshock' file_suffix])
    save([file_out '/Pss_tradeshock' file_suffix],'Pnow')
    
elseif nargin==2
    % file path supplied
    save([file_out '/ss_tradeshock'])
        save([file_out '/Pss_tradeshock'],'Pnow')
    
else
    % no file path supplied
    mkdir('mat_files')
    save('mat_files/ss_tradeshock')
    save('mat_files/Pss_tradeshock','Pnow')
end


