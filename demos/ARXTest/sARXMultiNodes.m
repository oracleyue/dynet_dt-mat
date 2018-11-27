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
simFlag  = 0;
plotFlag = 1;
% paths (use local drives instead of dropbox)
prePath  = [getenv('HOME') '/Workspace/swap/'];
dataPath = [prePath 'dynet/dataARX_p/'];
resPath  = [prePath 'dynet/resARX_p/'];
if ~exist(dataPath, 'dir') || ~exist(resPath, 'dir')
    mkdir(dataPath); mkdir(resPath);
end

% model
pList = [5 10 15 20];
m = 1;
SNR = 10;  % unit: dB
% dataset
nDataset = 50;
range_ts = 1:1000;
% solver
solver = 'gspmc';
% solver   = 'gsbl';     % sparsity method
lambda   = 1e-2;         % initial noise variance
% solver   = 'gil1';
% lambdaList = [5e-2 1e-1 1e-1 1e-1];
% model structure
norders  = [3 2 2];    % orders of ARX models
% file names
dataName = sprintf('arxsim_snr%undata%u', SNR,nDataset);

% ================================================================
%% Model Simulation
% ================================================================
if simFlag

    for p = pList
        fprintf('Model simulation: START\n')
        fprintf('    info: p%um%u SNR%u\n', p,m,SNR);
        simTimer = tic;

        % configure options
        [dataOpt, netOpt, signalOpt] = arxsimOptions(p,m);
        dataOpt.numData = 1000;
        dataOpt.numExpr = nDataset;        % #datasets
        dataOpt.idxStart = 200;            % #points discarded
        dataFName = [dataPath dataName ...
                      sprintf('p%um%u.mat', p,m)];
        dataOpt.fileName = dataFName;
        netOpt.sparsity = .2;              % sparsity density
        signalOpt.SNR = SNR;               % SNR in dB
        signalOpt.inputType = 'gauss';

        % simulation
        arx_sim(dataOpt, netOpt, signalOpt, 1);

        simEtime = toc(simTimer);
        fprintf('    elapsed time: %f sec\n', simEtime);
        fprintf('Model simulation: END.\n\n')
    end

end % END of simFlag

% ================================================================
%% Network Inference
% ================================================================
nP = length(pList);
% data struct
dynet_cells = cell(nDataset, nP);
perf_cells = cell(nDataset, nP);

for idxP = 1:nP
    p = pList(idxP);
    fprintf('Network reconstruction: START\n')
    fprintf('    info: p%um%u SNR%u\n', p,m,SNR);

    % load data
    dataFName = [dataPath dataName sprintf('p%um%u.mat', p,m)];
    load(dataFName);

    % extract model and sim info
    numExpr = nDataset;  % size(outputs,1);
    Ts = simInfo.Ts;
    num_ts = size(outputs{1},1);
    p = size(outputs{1},2);
    m = size(inputs{1},2);

    % process each dataset
    for idxExpr = 1:numExpr
        fprintf('    ... process Model %u \n', idxExpr);
        netkTimer = tic;

        % format data
        outdata = outputs{idxExpr}(range_ts,:);
        indata = inputs{idxExpr}(range_ts,:);
        data = iddata(outdata, indata, Ts);

        % choose labmda when using groupIRL1
        if strcmp(solver, 'gil1')
            lambda = lambdaList(idxP);
        end

        % network object
        netobj = idnet(p,m);
        solverOptions = solverInit(solver, lambda);
        netobj.arx(data, norders, solverOptions);

        % save ground truth in netobj
        netobj.GroundTruth.Q = models.Qb{idxExpr};

        netkEtime = toc(netkTimer);
        fprintf('        elapsed time: %f sec\n', netkEtime);

        % print results
        Qgt = netobj.GroundTruth.Q;
        Qq = netobj.BoolNet.Q;
        TP = sum(sum(Qgt & Qq));
        FP = sum(sum(~Qgt & Qq));
        FN = sum(sum(Qgt & ~Qq));
        TN = sum(sum(~Qgt & ~Qq));
        Prec = TP/(TP+FP);
        TPR  = TP/(TP+FN);
        fprintf('        Prec: %f%% \n', Prec*100);
        fprintf('        TPR: %f%% \n', TPR*100);

        % save results
        dynet_cells{idxExpr,idxP} = netobj;
        perf_cells{idxExpr,idxP} = [TP FP FN TN Prec TPR];
    end % END of numExpr

    fprintf('Network reconstruction: END\n\n')
end % END of pList
matName = [resPath sprintf('dynet_snr%undata%u_multiNodes.mat', ...
                           SNR,nDataset)];
save(matName, 'dynet_cells', 'perf_cells', 'pList', 'SNR');

% ================================================================
%% Visualization and Exportation
% ================================================================
if plotFlag

for idxP = 1:nP
    p = pList(idxP);
    for idxExpr = 1:numExpr
        netobj = dynet_cells{idxExpr, idxP};

        fighl = figure('visible', 'off');
        pg = plot(netobj);
        title(sprintf('Digraphs of Model %u (p=%u; SNR=%u dB)', ...
                      idxExpr, p, SNR));
        % highlight wrong links in red
        Qgt = logical(netobj.GroundTruth.Q);
        Qq = logical(netobj.BoolNet.Q);
        BoolNetWrong = ~Qgt & Qq;
        [s, t] = find(BoolNetWrong');
        highlight(pg, s', t', 'Edgecolor', 'r');
        % export figures in pdf
        figName = [resPath, sprintf('dynet_model%u_p%u.pdf', ...
                                    idxExpr, p)];
        set(fighl,'Units','Inches');
        pos = get(fighl,'Position');
        set(fighl,'PaperPositionMode','Auto', 'PaperSize',[pos(3), pos(4)]);
        print(fighl, figName, '-dpdf', '-r0')
    end
end

end % END of plotFlag
