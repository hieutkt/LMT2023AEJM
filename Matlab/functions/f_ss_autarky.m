% f_ss_autarky.m
%
%   This is the basic function to solve the model under autarky.
%
%   This function solves the model by varying P to clear C. There are three
%   inputs:
%       1. Input parameters
%       2. Location to export mat files (optional)
%       3. Suffix to tag onto file name (optional)

function output_moment = f_ss_autarky(par_in,file_out,file_suffix)

%% parameters and grids: set baseline parameter

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

%%% set initial grid size and prices
try 
    par.kstar  = par_in.new_kmax;
    Pnow        = par_in.new_Pnow;
catch
    disp('no new kmax /wstar provided');
    Pnow = 0.6530;
%     par2.kstar = 300;
%     Pnow = 0.8386;
    par.kstar = 12.7173;
end

%% Solve GE
   
% form par2 structure
par2 = par;

%%% Solve GE
P_guess = zeros(1000,1);
new_C_store = zeros(1000,1);
itermax = 1000;
tolw = 1e-3;
relx = 0.5;
for iter=1:itermax
    
    % set Q to P for now
    par2.P  = Pnow;
    if par2.lam==0 && par2.zet==0
        par2.Q  = par_in.Q;
    else
        par2.Q  = Pnow;
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
    % actual consumption
    new_Cd = (gr.y_i(:).^((par2.epsi-1)/par2.epsi)'*Px.g_ss(:))^(par2.epsi/(par2.epsi-1));
    new_C = new_Cd;
    % update value function
    initval.Ev = Ev;

    % check tolerance
    diff_C = new_C - old_C;
    if abs(diff_C) < tolw
        disp(diff_C);
        break
    end
    
    % update prices
    if iter==1
        
        new_C = min(100,new_C);
        
        lc = min(new_C,old_C);
        rc = max(new_C,old_C);
        
    else
        
        lc_cand = min(new_C,old_C);
        rc_cand = max(new_C,old_C);
        
        lc = max(lc,lc_cand);
        rc = min(rc,rc_cand);
        
    end
    if abs(lc-rc)< min(tolw,1e-6)
        lc = lc*0.999;
        rc = rc*1.001;
    end
    P_guess(iter) = 1/(((1-relx)*lc + relx*rc)*par2.chi);
%     P_guess(iter) = 1/(((1-relx)*old_C + relx*new_C)*par2.chi);
    Pnow = P_guess(iter);
    new_C_store(iter) = new_C;

    % plot progress
    subplot(311)
    drawnow
    hold on
    plot(iter,diff_C,'x');
    subplot(312)
    plot(iter,Pnow,'x');
    drawnow
    hold on
    subplot(313)
    plot(iter,new_C,'x');
    drawnow
    hold on

    % check grid sizes
    if par2.lam==0 && par2.zet==0
        % no need to check grid sizes if special case
        
%         disp(diff_C);
        
        if 1<abs(diff_C) && abs(diff_C)<10
            relx=.75;
        elseif 0.5<abs(diff_C) && abs(diff_C)<1
            relx = .5;
        elseif abs(diff_C)<0.5
            relx = .5;
        end
        
    else
        if max(pol.kpr_up(:)) < .9*gr.k_grid(end)
            % too big, truncate a little

            par2.kstar = .9*gr.k_grid(end);
            disp('new kmax')
            disp(num2str(par2.kstar))

        elseif max(pol.kpr_up(:)) >= gr.k_grid(end)
            % too small, expand a little

            par2.kstar = 1.1*gr.k_grid(end);
            disp('new kmax')
            disp(num2str(par2.kstar))

        else
            disp(diff_C);

            if 1<abs(diff_C) && abs(diff_C)<2
                relx=.75;
            elseif 0.5<abs(diff_C) && abs(diff_C)<1
                relx = .5;
            elseif abs(diff_C)<0.5
                relx = .25;
            end

        end
    end

end
close all

% compute moments
output_moment = get_moments(par2,N,gr,Px,pol,1);

% send these variables into the workspace to be saved and used in version
% where we "close" the economy
initval.Ev = Ev;
Pf = Pnow;
dM = 0;

% export full workspace
try
    if par_in.save == 1
        
        % save kmax for trade steady-state
        kmax_autarky = gr.k_grid(end);
        
        if nargin==3
            % file path and suffix supplied
            if par2.lam>0
                save([file_out '/ss_autarky' file_suffix])
                save([file_out '/Pss_autarky' file_suffix],'Pnow')
                save([file_out '/kmax_autarky' file_suffix],'kmax_autarky');
            else
                save([file_out '/ss_autarky_q1' file_suffix])
                save([file_out '/Pss_autarky_q1' file_suffix],'Pnow') 
                save([file_out '/kmax_autarky_q1' file_suffix],'kmax_autarky');
            end

        elseif nargin==2
            % file path supplied
            if par2.lam>0
                save([file_out '/ss_autarky'])
                save([file_out '/Pss_autarky'],'Pnow')
                save([file_out '/kmax_autarky'],'kmax_autarky');
            else
                save([file_out '/ss_autarky_q1'])
                save([file_out '/Pss_autarky_q1'],'Pnow') 
                save([file_out '/kmax_autarky_q1'],'kmax_autarky');
            end

        else
            % no file path supplied
            mkdir('mat_files')
            if par2.lam>0
                save('mat_files/ss_autarky')
                save('mat_files/Pss_autarky','Pnow')
                save('mat_files/kmax_autarky','kmax_autarky');
            else
                save('mat_files/ss_autarky_q1')
                save('mat_files/Pss_autarky_q1','Pnow') 
                save('mat_files/kmax_autarky_q1','kmax_autarky');
            end
        end
    end
catch
end




