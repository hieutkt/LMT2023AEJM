function bet = pop_reg(y,X,g,Px,N)
% pop_reg.m
%   This is a generic function to run population regressions
% y is the regressand of interest
% X are the regressors of interest
% g is the numerical approximation to the distribution
% Px, if supplied, is the conditional: g(t+1) = P(t+1|t)*g(t)

[xrows,xcols] = size(X);

% transpose into variables along the "row" dimension
if xrows>xcols
    X = X';
end
y= y(:);
g = g(:);
nstates = length(g);

% set rhs regressors
XX = [ones(1,nstates) ; X];

% attach probabilities to grid points
p_XX = bsxfun(@times,XX,g');

% E[xx']
EXX = XX*p_XX';
% E[xy]
if nargin==3
    %%% doing contemporaneous regressions
    
    % simply compute E[xy] normally
    E_XY = p_XX*y;
else
    %%% doing lagged regressions
    
    % first pull out current distribution
    g0 = Px.g_ss(:)/sum(Px.g_ss(:));
    p_XX_stay = bsxfun(@times,XX,g0');
    [nrows,~] = size(XX);
    pL_XX = zeros(size(XX));
    % iterate forward using conditionals
    for ii = 1:nrows
        xin = p_XX_stay(ii,:);
        temp = Px.H_s*reshape(Px.G_k*xin(:),N.s,N.k);
        pL_XX(ii,:) = temp(:);
    end

    % compute E[x(t)y(t+1)], and rescale if iterating leads to exits
    g1 = Px.H_s*reshape(Px.G_k*g0,N.s,N.k);
    m1 = sum(g1(:));
    E_XY = pL_XX*y/m1;
end

% get beta
bet = EXX\E_XY;


end