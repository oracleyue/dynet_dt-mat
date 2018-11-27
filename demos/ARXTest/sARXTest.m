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
simFlag = 0;
plotFlag = 1;
% paths
dataPath = './';
resPath  = './';

% model
p        = 10;         % num of outputs
m        = 1;          % num of inputs
SNR      = 10;
N        = 1000;       % length of ts in sim
% dataset
nDataset = 1;
range_ts = 1:1000;     % range of ts in use
% sparsity solver
% solver = 'gil1';     % lambda = .01
% solver = 'gsbl';
solver   = 'gspmc';
lambda   = .01;        % regul para or init noise variance
solverOptions = solverInit(solver, lambda);
% model structure
norders  = [3 2 2];    % orders (+1) of ARX models

% ================================================================
%% Model Simulation
% ================================================================
if simFlag

fprintf('Model simulation: START\n')
simTimer = tic;

% configure options
[dataOpt, netOpt, signalOpt] = arxsimOptions(p,m);
dataOpt.numData = N;
dataOpt.numExpr = nDataset;        % #datasets
dataOpt.idxStart = 200;            % #points discarded
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
    outdata = outputs{idxExpr}(range_ts,:);
    indata = inputs{idxExpr}(range_ts,:);
    data = iddata(outdata, indata, Ts);

    % network object
    netobj = idnet(p,m);
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
    fprintf('       Prec: %f%% \n', TP/(TP+FP)*100);
    fprintf('       TPR: %f%% \n', TP/(TP+FN)*100);
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
    filename = [figPath, 'dynet_model' num2str(idxExpr), '.pdf'];
    set(fighl,'Units','Inches');
    pos = get(fighl,'Position');
    set(fighl,'PaperPositionMode','Auto', 'PaperSize',[pos(3), pos(4)]);
    print(fighl, filename, '-dpdf', '-r0')
end

end
