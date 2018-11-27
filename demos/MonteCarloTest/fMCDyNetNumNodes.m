% function fMCDyNetNumNodes(SetIndex)
%
% FMCDYNETNUMNODES Run Monte Carlo study on different methods (GL:
% group Lasso; IRL: iterative reweighted l1 method; SBL: sparse
% Bayesian learning). In each run, 50 random networks are inferred
% using the specified methods.
%
% Inputs:
%   - SetIndex: integer
%     To specify #nodes of the random networks. In codes, the
%     specific datasets will be chosen.
%
% Copyright [2017] <oracleyue>



%% Environment settings
% path settings
addpath('./')
addpath('./auxiliary/');
addpath('./simulators/');
datpath = './SimData/';
respath = './Results/';

% control flags
sim_flag = 0;
inf_flag = 1;

% settings for solver
global glb_optimizer
glb_optimizer = 'cvx';
global MAX_ITER ABSTOL  %precision for prox
MAX_ITER = 1e4;
ABSTOL = 1e-8;
% algorithm settings
solver = 'gil1';
% lambda = .0005;
lambda = .00055;

% Monte Carlo settings
mc_rounds = 50;

% output variables
MC_Results = zeros(mc_rounds, 4);
  %  = [TypeI_Error, TypeII_Error, Correct_Percent, lambda]

% Start parpool
% parpool('local', 20);


%% Simluation settings
sys_info = struct;
sys_info.n = 80;
sys_info.p = 40;
numHiddenNodes = 1;
sys_info.IdxNodes = 1:(numHiddenNodes+1):sys_info.n;
sys_info.m = 1;
sys_info.N = 2000;
sys_info.TsFactor = 20*2;
sys_info.SNR = .0001;
sys_info.InitVal = 0;
sys_info.sparse_density = .005/4;


%% System simulation
if sim_flag
    sim_timer = tic;

    for ip = 1:mc_rounds
        series_prefix = ['NO.' sprintf('%03d',ip)];
        dynet_sim(sys_info, datpath, series_prefix);
    end

    sim_etime = toc(sim_timer);
    fprintf('The elapsed time for system simulation: %f seconds \n', ...
        sim_etime);
end


%% Network inference
if inf_flag
net_timer = tic;

for ip = 1:mc_rounds

% start timer
netk_timer = tic;

% % Load datasets
% path and names
series_prefix = ['NO.' sprintf('%03d',ip)];
fname_prefix = ['dynetsim_' series_prefix '_' num2str(sys_info.p) ...
    '_' num2str(sys_info.SNR)];
% dimensions
p = sys_info.p;
m = sys_info.m;
SNR = sys_info.SNR;
% reading from txt files
fname = [datpath fname_prefix '_out.txt'];
output = dlmread(fname, '\t');
fname = [datpath fname_prefix '_in.txt'];
input = dlmread(fname, '\t');
fname = [datpath fname_prefix '_info.txt'];
sim_info = dlmread(fname, '\t');
Qgt = sim_info(1:p, 1:p);
Ts_base = sim_info(p+1,1);
clear sim_info


% % Show the ground truth of digraphs
% % plot the ground truth
% fig_hl = figure('visible', 'off');
% set(fig_hl, 'units', 'inches', 'position', [18.4444 7.8333 7.7778 5.8333]);
% Ggt = digraph(Qgt');
% pg = plot(Ggt, 'Layout', 'force');
% pg.NodeColor = 'red';
% title('Ground truth of Boolean network')

% % Distill connected components
% uG = graph(Qgt | Qgt');
% uG = digraph(Qgt);
% bins = conncomp(uG);
% for i = unique(bins)
%     countGroup(i) = length(find(bins==i));
% end
% [val, idx] = max(countGroup);
% nodeSet = find(bins == idx);


% Setup datasets
%   Notes: DT = Ts_sysali/10 /Ts_base; DT = floor(DT);
% DT = 2;
DTList = [2, 4, 6, 8, 12]; DT = DTList(SetIndex);
% NsamplesList = [300 500 700 900 1100]; Nsamples = NsamplesList(SetIndex);
Nsamples = 1100;
Ts = Ts_base * DT;
output = output(101:DT:Nsamples, :);
input = input(101:DT:Nsamples, :);

% % Select nodes from a connected component in graph:
% select_index = 1:p;
% % select_index = nodeSet;
% output = output(:, select_index);
% input = input(:, m); % input = zeros(size(output,1), 1);
% p = length(select_index);
% Qgt = Qgt(select_index, select_index);

% Building data structures
data = iddata(output, input, Ts);
netobj = idnet(p,m);
norders = {ones(p,1)*2, ones(p,p)*2, ones(p,m)*1, []};
clear input output

% Network inference
netobj.ThresholdZero = 1e-6; % used for generating bool network
% solver = 'gil1';
Lambda = [0, lambda, 0]; %gil1:.001; gl1:2
netobj.arx(data, norders, {solver, Lambda});

clear data norders

% boolnet = netobj.bool(1e-6);
% Print results
% dynet = netobj.dsf();


%% Publish digraphs in pdf
fig_hl = figure('visible', 'off');
pg = plot(netobj);
xlabeltext = sprintf(['The sampling period = ' num2str(DT) 'Ts (second);\n ' ...
                    'The samples: 100:' num2str(DT) ':' num2str(Nsamples) '; '...
                    'The noise strength: ' num2str(SNR) ';\n'...
                    'The threshold of zero in Boolean network = ' num2str(netobj.ThresholdZero) '.']);
xlabel(xlabeltext)
title([SolverNameMap(solver) ' (\lambda = ' num2str(lambda) ')']);
% highlight wrong links in red
Qq = logical(netobj.BoolStruct.Q);
BoolNetWrong = ~Qgt & Qq;
[s, t] = find(BoolNetWrong');
highlight(pg, s', t', 'Edgecolor', 'r');

% publish figures in pdf
figpath = [respath 'outputs/'];
figname_prefix = strrep(fname_prefix, 'dynetsim', 'dynet');
filename = [figpath, figname_prefix, '_' solver, '.pdf'];
set(fig_hl,'Units','Inches');
pos = get(fig_hl,'Position');
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0')

clear netobj

%% Publish Monte Carlo results in pdf
num_all_links_gt = sum(sum(Qgt));
num_all_links_res = sum(sum(Qq));
num_correct_links = sum(sum(Qgt & Qq));
num_wrong_links = sum(sum(~Qgt & Qq));
num_missed_links = sum(sum(Qgt & ~Qq));

MC_Results(ip,4) = lambda;
MC_Results(ip,1) = num_wrong_links / num_all_links_res;
MC_Results(ip,2) = num_missed_links / num_all_links_gt;
MC_Results(ip,3) = num_correct_links / num_all_links_gt;


% time benchmark
netk_etime = toc(netk_timer);
fprintf('The elapsed time of inferce No.%d: %f seconds \n', ...
        ip, netk_etime);
fprintf('Type I  Error: %f%% \n', MC_Results(ip,1)*100);
fprintf('Type II Error: %f%% \n\n', MC_Results(ip,2)*100);

end


%% saving results in MAT
fname_prefix = ['MC' num2str(mc_rounds) '_p' ...
                num2str(sys_info.p) '_SNR' num2str(sys_info.SNR)...
                '_Ts' num2str(DT)...
                '_Ns' num2str(Nsamples)];
txtname = [respath fname_prefix '.txt'];
fid = fopen(txtname, 'w');
fspec = '%f \t %f \t %f \t %f\n';
fprintf(fid, fspec, MC_Results');
fclose(fid);

% time benchmark
net_etime = toc(net_timer);
fprintf('The elapsed time of network inferce: %f seconds \n', ...
        net_etime);

% shutdown parpool
% delete(gcp);

end   %END: if inf_flag

quit
