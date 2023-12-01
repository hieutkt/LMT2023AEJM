clc
clear

in_path = 'mat_files_compiled/baseline/mat_files/results/';

mkdir('../stata_replication/Data/output/')
out_path_stata = '../stata_replication/Data/output/';

load([in_path 'stats_data.mat'])

ind = [1:32]';
year = [1:15]';

A =ones(32,1);
B =ones(15,1);

ind2 = kron(ind,B);
year2 = kron(A,year);

tfp = reshape(stats_agg.TFPQ_t,[15*32,1]);
sd_mpkt = reshape(stats_agg.sd_mpk_t,[15*32,1]);
inact = reshape(stats_agg.inact_t,[15*32,1]);
ipos = reshape(stats_agg.ipos_t,[15*32,1]);
ineg = reshape(stats_agg.ineg_t,[15*32,1]);
avg_ik = reshape(stats_agg.avg_ik_t,[15*32,1]);
agg_ik = reshape(stats_agg.agg_ik_t,[15*32,1]);
exit_rate = reshape(stats_agg.exit_rate_t,[15*32,1]);

imppen = reshape(stats_agg.import_pen_diff_t,[15*32,1]);
imppen_r = reshape(stats_agg.import_pen_t_realized,[15*32,1]);

DATA = [ind2 year2 tfp sd_mpkt inact ipos ineg avg_ik agg_ik exit_rate imppen imppen_r];

writematrix(DATA,[out_path_stata 'dataforstata_agg.xlsx'])

%-----------------------------------------------
%microdata - inaction
clc
load([in_path 'stats_data.mat'])

ind2 = stats_micro.micro1.ind_tag;
year2 = stats_micro.micro1.year_tag;
id = stats_micro.micro1.id_tag;

logk = stats_micro.micro1.logk_data;
logs = stats_micro.micro1.logs_data;
inv_rate = stats_micro.micro1.inv_rate_data;
inact = stats_micro.micro1.inact_data;
logmrpk = stats_micro.micro1.log_mpk_data;

imppen = stats_micro.micro1.imp_diff_data;
imppen_r = stats_micro.micro1.imp_pen_data;

DATA2 = [ind2 year2 id logk logs inv_rate inact logmrpk imppen imppen_r];

writematrix(DATA2,[out_path_stata 'dataforstata_micro1.csv'])

%-----------------------------------------------

%microdata - exit
clc
load([in_path 'stats_data.mat'])

ind2 = stats_micro.micro2.ind_tag;
year2 = stats_micro.micro2.year_tag;
id = stats_micro.micro2.id_tag;

logk = stats_micro.micro2.logk_data;
logs = stats_micro.micro2.logs_data;
stay = stats_micro.micro2.i_stay_data;

imppen = stats_micro.micro2.imp_diff_data;
imppen_r = stats_micro.micro2.imp_pen_data_2;


DATA2 = [ind2 year2 id logk logs stay imppen imppen_r];

writematrix(DATA2,[out_path_stata 'dataforstata_micro2.csv'])

%-----------------------------------------------

%frictionless
clc
load([in_path 'stats_frictionless_data.mat'])

ind2 = stats_micro.micro1.ind_tag;
year2 = stats_micro.micro1.year_tag;
id = stats_micro.micro1.id_tag;

logk = stats_micro.micro1.logk_data;
logs = stats_micro.micro1.logs_data;
inv_rate = stats_micro.micro1.inv_rate_data;
inact = stats_micro.micro1.inact_data;
logmrpk = stats_micro.micro1.log_mpk_data;

imppen = stats_micro.micro1.imp_diff_data;
imppen_r = stats_micro.micro1.imp_pen_data;

DATA2 = [ind2 year2 id logk logs inv_rate inact logmrpk imppen imppen_r];

writematrix(DATA2,[out_path_stata 'dataforstata_frict.csv'])

%-----------------------------------------------

%microdata - exit
clc
load([in_path 'stats_frictionless_data.mat'])

ind2 = stats_micro.micro2.ind_tag;
year2 = stats_micro.micro2.year_tag;
id = stats_micro.micro2.id_tag;

logk = stats_micro.micro2.logk_data;
logs = stats_micro.micro2.logs_data;
stay = stats_micro.micro2.i_stay_data;

imppen = stats_micro.micro2.imp_diff_data;
imppen_r = stats_micro.micro2.imp_pen_data_2;

DATA2 = [ind2 year2 id logk logs stay imppen imppen_r];

writematrix(DATA2,[out_path_stata 'dataforstata_micro2_frict.csv'])

%-----------------------------------------------

clearvars -except in_path out_path_stata

load([in_path 'stats_frictionless_data.mat'])

ind = [1:32]';
year = [1:15]';

A =ones(32,1);
B =ones(15,1);

ind2 = kron(ind,B);
year2 = kron(A,year);

tfp = reshape(stats_agg.TFPQ_t,[15*32,1]);
sd_mpkt = reshape(stats_agg.sd_mpk_t,[15*32,1]);
inact = reshape(stats_agg.inact_t,[15*32,1]);
ipos = reshape(stats_agg.ipos_t,[15*32,1]);
ineg = reshape(stats_agg.ineg_t,[15*32,1]);
avg_ik = reshape(stats_agg.avg_ik_t,[15*32,1]);
agg_ik = reshape(stats_agg.agg_ik_t,[15*32,1]);
exit_rate = reshape(stats_agg.exit_rate_t,[15*32,1]);

imppen = reshape(stats_agg.import_pen_diff_t,[15*32,1]);
imppen_r = reshape(stats_agg.import_pen_t_realized,[15*32,1]);

DATA = [ind2 year2 tfp sd_mpkt inact ipos ineg avg_ik agg_ik exit_rate imppen imppen_r];

writematrix(DATA,[out_path_stata 'dataforstata_agg_fric.xlsx'])
