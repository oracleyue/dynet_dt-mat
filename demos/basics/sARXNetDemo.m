% This script generates and simulates sparse stable ARX models, and
% perform network inference.

% Copyright [2018] <oracleyue>
% Last modified on 19 Jul 2018


clear all
close all

% ================================================================
%% Environment Setup
% ================================================================

% setup search path
run('../../dynet_startup.m')
% flag
simFlag = 1;
plotFlag = 1;
% paths
dataPath = './';
resPath = '../../Results/resARX/';

% model
p = 10;
m = 1;
SNR = 10;
% dataset
nDataset = 1;
% solver
solver   = 'gsbl';     % sparsity method
lambda   = 1e-2;       % regularization parameter
norders  = [3 2 2];    % orders of ARX models


% ================================================================
%% Model Simulation
% ================================================================
if simFlag

fprintf('Model simulation: START\n')
simTimer = tic;

% configure options
[dataOpt, netOpt, signalOpt] = arxsimOptions();
dataOpt.numData = 1000;
dataOpt.numExpr = nDataset;        % #datasets
dataOpt.idxStart = 200;            % #points discarded
netOpt.numNodes = p;               % dim of network
netOpt.sparsity = .2;              % sparsity density
signalOpt.SNR = SNR;               % SNR in dB
signalOpt.inputType = 'gauss';

% simulation
arx_sim(dataOpt, netOpt, signalOpt, 1);

simEtime = toc(simTimer);
fprintf('    elapsed time: %f sec\n', simEtime);
fprintf('Model simulation: END.\n\n')

end

% ================================================================
%% Network Inference
% ================================================================
fprintf('Network reconstruction: START\n')

% load data
load([dataPath 'ARXSimData.mat'])

% extract model and sim info
numExpr = size(outputs,1);
Ts = simInfo.Ts;
num_ts = size(outputs{1},1);
p = size(outputs{1},2);
m = size(inputs{1},2);

% data struct
dynet_cells = cell(numExpr,1);

% process each dataset
for idxExpr = 1:numExpr
    fprintf('   ... process Model %u \n', idxExpr);
    netkTimer = tic;

    % format data
    data = iddata(outputs{idxExpr}, inputs{idxExpr}, Ts);

    % network object
    netobj = idnet(p,m);
    solverOptions = solverInit(solver, lambda);
    netobj.ThresholdZero = 1e-6;
    netobj.arx(data, norders, solverOptions);

    % save ground truth in netobj
    netobj.GroundTruth.Q = models.Qb{idxExpr};

    % save dynet
    dynet_cells{idxExpr} = netobj;

    netkEtime = toc(netkTimer);
    fprintf('       elapsed time: %f sec\n', netkEtime);

    % print results
    Qgt = netobj.GroundTruth.Q;
    Qq = netobj.BoolNet.Q;
    nLinkGT = sum(sum(Qgt));
    nLinkRes = sum(sum(Qq));
    TP = sum(sum(Qgt & Qq));
    FP = sum(sum(~Qgt & Qq));
    FN = sum(sum(Qgt & ~Qq));
    TN = sum(sum(~Qgt & ~Qq));
    fprintf('       Prec: %f%% \n', TP/(TP+FP));
    fprintf('       TPR: %f%% \n', TP/(TP+FN));
end

fprintf('Network reconstruction: END\n\n')

% ================================================================
%% Visualization and Exportation
% ================================================================
if plotFlag

for idxExpr = 1:numExpr
    netobj = dynet_cells{idxExpr};

    fighl = figure;
    pg = plot(netobj);
    title(sprintf('Digraphs of Model %u', idxExpr));
    % highlight wrong links in red
    Qgt = logical(netobj.GroundTruth.Q);
    Qq = logical(netobj.BoolNet.Q);
    BoolNetWrong = ~Qgt & Qq;
    [s, t] = find(BoolNetWrong');
    highlight(pg, s', t', 'Edgecolor', 'r');
    % export figures in pdf
    figPath = [resPath ''];
    % filename = [figPath, fignamePrefix, '_' solver, '.pdf'];
    filename = [figPath, 'dynet_model_' num2str(idxExpr), '.pdf'];
    pos = get(fighl,'Position');
    set(fighl,'PaperPositionMode','Auto', 'PaperSize',[pos(3), pos(4)]);
    print(fighl, filename, '-dpdf', '-r0')
end

end
