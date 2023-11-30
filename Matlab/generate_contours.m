% generate_contours.m
%
%   Generates contour plots for paper

clc
clear
close all

% data source folder
source_path = ['../stata_replication/Data/Modified/'];

% figure folder
figpath = 'figures/';

% four variants
suf = {'','_ip','_ipiv','_iv'};

% redefine figure colors for consistency (new matlab versions have dif
% color)
xblue = [0, 0.4470, 0.7410];
xred = [0.8500, 0.3250, 0.0980];

%% New hybrid data (May 2020)

kmin    = exp(9.19);
kmax    = exp(18.1);
smin    = 5.7;
smax    = 9.7;
ix      = 4;

%%% read in margins generated from stat
dat_in = xlsread([source_path 'hybrid_margins_data' suf{ix} '.xlsx']);

% base grid (from data)
log_cap_type_tot0 = dat_in(:,4);
tfp_va_k0 = dat_in(:,3);
try
    tfpv = linspace(smin,smax,100);
    capv = linspace(log(kmin),log(kmax),100);
catch
    tfpv = linspace(min(tfp_va_k0),max(tfp_va_k0),100);
    capv = linspace(log(10000),max(log_cap_type_tot0),100);
end

% find tfp consistent with pr
np = 1;
%     pr_fixed = linspace(0.40,0.9,np);
pr_fixed = 0.5;

% form interpolant for pr0
survive0 = dat_in(:,2);
Fs0 = scatteredInterpolant(tfp_va_k0,log_cap_type_tot0,survive0);
    % solve for pr0
tfpv0_star = zeros(100,np);
for ipr = 1:length(pr_fixed)
    diff_pr = @(tfp_in) -1*(Fs0(tfp_in,capv(:)) - pr_fixed(ipr)).^2;
    tfpv0_star(:,ipr) = goldenx(diff_pr,tfp_va_k0(1)*ones(size(capv(:))),tfp_va_k0(end)*ones(size(capv(:))));
end

% form interpolant for pr1
survive1 = dat_in(:,5);
Fs1 = scatteredInterpolant(tfp_va_k0,log_cap_type_tot0,survive1);
    % solve for pr1
tfpv1_star = zeros(100,np);
for ipr = 1:length(pr_fixed)
    diff_pr = @(tfp_in) -1*(Fs1(tfp_in,capv(:)) - pr_fixed(ipr)).^2;
    tfpv1_star(:,ipr) = goldenx(diff_pr,tfp_va_k0(1)*ones(size(capv(:))),tfp_va_k0(end)*ones(size(capv(:))));
end

%%% Plot outputs

% baseline and predicted survivial probs
figure
hold on
h1 = plot(capv,tfpv0_star,'Color',xblue,'linewidth',3);
h2 = plot(capv,tfpv1_star,'--','Color',xblue,'linewidth',3);
hold off
ylabel('log $\omega$','interpreter','latex','fontsize',24)
xlabel('log Capital Stock','interpreter','latex','fontsize',24)
grid on
    % user selects axis (auto enters after 5s)
t = timer('ExecutionMode', 'singleShot', 'StartDelay', 5, 'TimerFcn', @pressEnter);
start(t);
input_key = input('Input axis as [xmin xmax ymin ymax]: ');
stop(t);
delete(t);
if isempty(input_key)
    disp('Default axis');
else
    axis(input_key);
end
    % user selects legend location
t = timer('ExecutionMode', 'singleShot', 'StartDelay', 5, 'TimerFcn', @pressEnter);
start(t);
input_key2 = input('Input legend location (northwest, etc): ','s');
stop(t);
delete(t);
if isempty(input_key2)
    disp('Default legend (southwest)');
    input_key2 = 'southwest';
end
legend({'No shock','1 s.d. shock'},'fontsize',24,'location',input_key2)
print([figpath '/s_probs_china_shock_contour' suf{ix}],'-depsc')


%%
function pressEnter(HObj, event)
  import java.awt.*;
  import java.awt.event.*;
  rob = Robot;
  rob.keyPress(KeyEvent.VK_ENTER)
  rob.keyRelease(KeyEvent.VK_ENTER)
end
