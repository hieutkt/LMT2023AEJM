function ii = sim_discrete_iid(PP,PS,type)
% sim_discrete_iid.m. sim_discrete_iid simulates IID draws.
%
%   ii = sim_discrete_iid[PP,PS) generates the index for a univariate IID
%   process. PP either the pdf or the cdf (default). and PS is the number
%   of iid draws to be generated

% generate iid random uniform draws
UU = rand(PS,1);

if nargin==3
    switch type
        case('pdf')
            PP = cumsum(PP(:)');
        case('cdf')
            PP = PP(:)';
    end
end
        
% draw
ii = max(bsxfun(@times,bsxfun(@ge, UU(:), PP),1:length(PP)),[],2)+1;
