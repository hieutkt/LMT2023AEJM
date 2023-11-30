function k_grid = gen_k_grid(k_max,delta,del_grid,k_n)
% construct grid for k

lk_max = log(k_max);
 
lk_min = log(k_max)+(k_n-1)/del_grid*log(1-delta);
 
lk_grid = linspace(lk_min, lk_max , k_n)';
 
k_grid = exp(lk_grid);