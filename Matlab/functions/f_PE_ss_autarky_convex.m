% f_PE_ss_autarky_convex.m
%
%   This is the basic function to solve the model under autarky. 
%
%   This function solves the model by varying M to clear C. There are three
%   inputs:
%       1. Input parameters
%       2. Location to export mat files (optional)
%       3. Suffix to tag onto file name (optional)

function output_moment = f_PE_ss_autarky_convex(par_in,file_out,file_suffix)

%% parameters and grids: set baseline parameter

disp("[f_PE_ss_autarky_convex.m] Setup parameters and grids...")

% retrieve baseline parameters and grid sizes
output  = mainfuncMP_convex([],[],[],[],1,par_in);
par     = output.par;
N       = output.N;
clear output

% set initial max grid size and prices
par.kstar   = par_in.new_kmax;
Pnow        = par_in.new_Pnow;


%% Solve GE


disp("[f_PE_ss_autarky_convex.m] Solving individuals' problem...")
disp("    We're repeatedly solve until k grid is of appropriate size...")

% form par2 structure to avoid changing parameter input structure
par2 = par;

% Impose Q = P
par2.P  = Pnow;
par2.Q  = Pnow;

% Solve individual's problem
    % repeatedly solve until k grid is of appropriate size
i_start = 1;
i_too_big = 0;
i_too_small = 0;
n_count = 0;
while i_too_big || i_too_small || i_start
    
    % initialize
    i_start = 0;
    clear gr

    %%% set the grid
    gr = mainfuncMP_convex(par2,N,[],[],2);

    %%% solve the individual's problem
    try
        output  = mainfuncMP_convex(par2,N,gr,[],3,initval);
    catch
        output  = mainfuncMP_convex(par2,N,gr,[],3);
    end
    pol     = output.pol;
    Ev       = output.Ev;
    clear output

    % check grid sizes, but stop if iterating >50
    if n_count>50
        break
    end
    if max(pol.kpr_up(:)) < .9*gr.k_grid(end)
        % too big, truncate a little

        par2.kstar = .9*gr.k_grid(end);

        i_too_big = 1;
        i_too_small = 0;
        n_count = n_count + 1;

    elseif max(pol.kpr_up(:)) >= 0.999*gr.k_grid(end)
        % too small, expand a little

        par2.kstar = 1.1*gr.k_grid(end);

        i_too_big = 0;
        i_too_small = 1;
        n_count = n_count + 1;

    else
        % appropriate size reached
        
        i_too_big = 0;
        i_too_small = 0;

    end

end

% get distribution
output  = mainfuncMP_convex(par2,N,gr,pol,4);
Px      = output;
clear output

% Compute consumption given guess of P to get initial guess for solver
old_C = 1/(par2.P*par2.chi);   
new_Cd = (gr.y_i(:).^((par2.epsi-1)/par2.epsi)'*Px.g_ss(:))^(par2.epsi/(par2.epsi-1));
Me0 = (old_C/new_Cd)^((par2.epsi-1)/par2.epsi); % initial guess

% GE step
if par_in.cali == 1
    % if calibrating, moments of interest invariant to GE. skip this.
    Me_new = 1;
    
else
    % find scaling factor such that new_Cd = old_C
    Xscale = fsolve(@(x) Fdiff_C(x, par2, N, gr, pol, old_C),par2.Me*Me0);
    Me_new = Xscale*par2.Me;
    
end

%% Get statistics here: solve again with new scaling factor

% set the grid
par2.Me = Me_new;
gr = mainfuncMP_convex(par2,N,[],[],2);

% solve the individual's problem
try
    output  = mainfuncMP_convex(par2,N,gr,[],3,initval);
catch
    output  = mainfuncMP_convex(par2,N,gr,[],3);
end
pol     = output.pol;
Ev       = output.Ev;
clear output

% get distribution
output  = mainfuncMP_convex(par2,N,gr,pol,4);
Px      = output;
clear output
% get new consumption to check
new_Cd = (gr.y_i(:).^((par2.epsi-1)/par2.epsi)'*Px.g_ss(:))^(par2.epsi/(par2.epsi-1));

% compute moments
output_moment = get_moments_convex(par2,N,gr,Px,pol,1);
output_moment.avg.new_Cd = new_Cd;
output_moment.avg.old_C = old_C;
output_moment.par.Me_new = Me_new;

% send these variables into the workspace to be saved
initval.Ev = Ev;
Pf = Pnow;
dM = 0;
if par_in.save == 1
    
    % save kmax for trade steady-state
    kmax_autarky = gr.k_grid(end);
    
    % export full workspace
    if nargin==3
        % file path and suffix supplied
        save([file_out '/ss_autarky' file_suffix])
        save([file_out '/Pss_autarky' file_suffix],'Pnow')
        
        save([file_out '/kmax_autarky' file_suffix],'kmax_autarky');
        
    elseif nargin==2
        % file path supplied
        save([file_out '/ss_autarky'])
        save([file_out '/Pss_autarky'],'Pnow')
        
        save([file_out '/kmax_autarky'],'kmax_autarky');

    else
        % no file path supplied
        mkdir('mat_files')
        save('mat_files/ss_autarky')
        save('mat_files/Pss_autarky','Pnow')
        
        save('mat_files/kmax_autarky','kmax_autarky');
        
    end

end


end

%% Function to compute excess demand
function  yy = Fdiff_C(Mx,par, N, gr, pol, old_C)
    gr_new = gr;
    gr_new.ps_dist = gr.ps_dist*Mx;
    output_new  = mainfuncMP_convex(par,N,gr_new,pol,4);
    new_Cd = (gr.y_i(:).^((par.epsi-1)/par.epsi)'*output_new.g_ss(:))^(par.epsi/(par.epsi-1));
    yy = old_C - new_Cd;

end

% iter=0;
% for Mx = 0.1:.1:5
%     iter = iter+1;
%     gr_new = gr;
%     gr_new.ps_dist = gr.ps_dist*Mx;
%     output_new  = mainfuncMP_convex(par2,N,gr_new,pol,4);
%     new_Cd = (gr.y_i(:).^((par2.epsi-1)/par.epsi)'*output_new.g_ss(:))^(par2.epsi/(par2.epsi-1));
%     yy(iter) = old_C - new_Cd;
% end
