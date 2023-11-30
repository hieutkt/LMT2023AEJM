function [ YY ] = sim_discreteMC(PP,T,t_burn,PS,init_conds)
% sim_discreteMC.m
%
% Simulates a system of discrete Markov chain, given
%   1. A markov transition matrix PP
%   2. The # of variables in the simulation PS
%   3. Simulation length T
%   4. # of people in panel
%   5. Initial conditions
%
% Notes:
%   1. The simulation only returns the indices corresponding to the matrix PP;
%       user has to use the index to recover the actual value
%   2. The simulation returns the entire series along with the burn in
%   time. User has to remove burn in period outside of function

% pre-allocate memory
YY = zeros(T+t_burn+1,PS);

% generate iid random uniform draws
UU = rand(T+t_burn+1,PS);

% generate index
[~,ncols] = size(PP);
ind_vec = repmat(1:ncols,PS,1);

%== simulate forward markov chain
% initialize process
PP_cdf = cumsum(PP,2); % conditional CDF
[init_conds_II, col_vec_II] = ndgrid(init_conds,1:ncols);
PP_ind_II = sub2ind(size(PP_cdf),init_conds_II,col_vec_II);
ii = max(bsxfun(@ge, UU(1,:)', PP_cdf(PP_ind_II)).*ind_vec,[],2)+1; %+1 for correction
YY(1,:) = ii;
% simulate for T periods
for t=2:T+t_burn+1
    [ii_II, col_vec_II] = ndgrid(ii,1:ncols);
    PP_ind_II = sub2ind(size(PP_cdf),ii_II,col_vec_II);    
    ii = max(bsxfun(@ge, UU(t,:)', PP_cdf(PP_ind_II)).*ind_vec,[],2)+1;
    YY(t,:) = ii;
end

