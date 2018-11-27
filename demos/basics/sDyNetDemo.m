% Perform network inference on a sepecific dataset

% Copyright (c) 2014-2018, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 20 Jul 2018


close all
clear all

% ================================================================
%% Environment settings
% ================================================================

% path settings
run('../../dynet_startup.m')
dataPath = '../../SimData/';
resPath = '../../Results/';

% flags
simFlag  = 0;
compFlag = 0;

% network
p        = 10;         % dim of output
m        = 1;          % dim of input
SNR      = 60;         % signal-noise-ration
sigma    = 1;          % variance when Gaussian; magnitude when step
% datasets
nDataset = 1;          % num of datasets to generate
dataID   = 1;          % index of dataset used
% solver
solver   = 'gil1';     % sparsity method
lambda   = 3e-3;     % regularization parameter
norders  = [2 5 5];    % orders of ARX models


% ================================================================
%% Model Simulation
% ================================================================
if simFlag
simTimer = tic;

% model setup
model.p = p;                            % dim of outputs
model.n = model.p * 2;                  % dim of states
% - choose hidden states: fixed equal number in each edge
%             numHiddenNodes = 2;
%             model.IdxNodes = 1:(numHiddenNodes+1):model.n;
% - choose hidden states: randomly select
%             indexPerm = randperm(model.n);
%             model.IdxNodes = indexPerm(1:model.p);
% - choose hidden states: randomly select; >= 1 in each edge
indexPerm = randperm(model.n);
indexPermOddIndices = find(rem(indexPerm,2) ~= 0);
indexPermOdd = indexPerm(indexPermOddIndices);
model.IdxNodes = indexPermOdd(1:model.p);
model.m = m;                            % dim of inputs
model.N = 1000;                         % #samples
model.TsRatio = 80;                     % see doc of "dynet_sim.m"
model.sigma = sqrt(sigma^2 / 10^(SNR/10));
model.x0 = 0;                           % initial value of states
model.SparseDensity = .02/p;
% inputSignal.type = 'gauss';
% inputSignal.param = [0, sigma];
            inputSignal.type = 'step';
            inputSignal.param = 1;
pathName = '../../SimData/'; %  path to save data/figures
for indexJob = 1:nDataset
    fileName = ['dynet' '_p' num2str(model.p) '_SNR' num2str(SNR) 'dB'...
                '_Data' sprintf('%02d', indexJob) ];
    dynet_sim(model, 'PathName', pathName, 'FileName', fileName, ...
              'Input', inputSignal, 'Mode', 'measurement', ...
              'Export', 0);
end
% save data
fid = fopen([pathName 'dataList.log'], 'a');
fprintf(fid, [fileName '\n']);
fclose(fid);

simEtime = toc(simTimer);
fprintf('The elapsed time for system simulation: %f seconds \n', ...
        simEtime);
end


% ================================================================
%% Network Reconstruction
% ================================================================
% Start parpool
% parpool('local', 4);

netTimer = tic;

% load datasets
fnamePrefix = ['dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_Data'...
               sprintf('%02d',dataID)];
load([dataPath 'dataMAT/' fnamePrefix '_info.mat']);  % Qgt, Ts_sim
load([dataPath 'dataMAT/' fnamePrefix '.mat']);       % input, output, Ts

% data struct
data = iddata(output.y, input, Ts);

% reconstruction
netobj = idnet(p,m);
netobj.ThresholdZero = 1e-6;
solverOptions = solverInit(solver, lambda);
netobj.arx(data, norders, solverOptions);


% ================================================================
%% Visualization and Exportation
% ================================================================
% fig_hl = figure('visible', 'off');
fig_hl = figure;
pg = plot(netobj);
solverName = SolverNameMap(solver);
xlabel([solverName.abbrName ' (\lambda = ' num2str(lambda) ')']);
% highlight wrong links in red
Qq = logical(netobj.BoolNet.Q);
BoolNetWrong = ~Qgt & Qq;
[s, t] = find(BoolNetWrong');
highlight(pg, s', t', 'Edgecolor', 'r');

% publish figures in pdf
figPath = [resPath ''];
filename = [figPath, fnamePrefix, '_' solver, '_Lambda_' num2str(lambda), '.pdf'];
pos = get(fig_hl,'Position');
set(fig_hl,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0')

% print results
nLinkGT = sum(sum(Qgt));
nLinkRes = sum(sum(Qq));
TP = sum(sum(Qgt & Qq));
FP = sum(sum(~Qgt & Qq));
FN = sum(sum(Qgt & ~Qq));
TN = sum(sum(~Qgt & ~Qq));

netEtime = toc(netTimer);
fprintf('Elapsed time on dynet No.%d: %f seconds \n', ...
        dataID, netEtime);
fprintf('Precision: %f%% \n\n', TP/(TP+FP));
fprintf('TPR: %f%% \n', TP/(TP+FN));


% ================================================================
%% Model Comparison
% ================================================================

if compFlag

dsfobj = dynet(netobj.dsf());
[fitScores, yModelResponse, xInitCond] = ...
    dsfobj.compare(data, 'plot');

end   %END: if compFlag


% shutdown parpool
% delete(gcp);
