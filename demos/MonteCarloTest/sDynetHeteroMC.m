% Perform network inference on one hetergeneous dataset, which consists
% of input-output time series of multiple experiments.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 23 Aug 2017



%% Environment Settings

% % ==========================================================
% % NOTES: comment this section if use "batchproc.m"
% % ==========================================================
% % path settings
% run('../../dynet_startup.m')
% dataPath = '../../SimData/MCDataSNRs/dataMAT/';
% resPath  = '../../Results/';

% % solver settings
% % solverOptions = solverInit('gsbl', 9e-3, 0);
% solverOptions = solverInit('gil1', 5.5e-3);

% % dataset settings
% p          = 10;     % dim of y
% m          = 1;      % dim of u
% SNR        = 40;     % signal-to-noise ratio
% dataIDList = 1:5;    % dataset ID No
% % ==========================================================

% variables: memory, I/O, files
netResults = zeros(5, length(dataIDList));
solverName = SolverNameMap(solverOptions.solver);
txtName = [resPath 'dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_' ...
           solverName.abbrName '.txt'];
% if exist(txtName, 'file'), delete(txtName); end
fid = fopen(txtName, 'w');
fspec = 'TP\tFP\tFN\tTN\tlambda\n';  % print the column names for tables
fprintf(fid, fspec); fclose(fid);


%% Network Inference

% Start parpool and timer
% parpool('local', 4);
netTimer = tic;

for indexDataID = 1:length(dataIDList)
    dataID = dataIDList(indexDataID);

    % Load Datasets
    fnamePrefix = ['dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_Data'...
                   sprintf('%02d',dataID)];
    load([dataPath fnamePrefix '_info.mat']);  % var: Qgt, Ts_sim
    load([dataPath fnamePrefix '.mat']);       % var: input, cOutputs, Ts

    % Build iddata
    data = {};
    for k = 1:length(cOutputs)
        output = cOutputs{k};
        data{k} = iddata(output, input, Ts);
    end

    % Parameter estimations
    netobj = idnet(p,m);
    norders = [2 5 5];
    netobj.ThresholdZero = 1e-6; % used for generating bool network
    netobj.arx(data, norders, solverOptions);
    boolnet = netobj.BoolNet;    % netobj.bool('L', 1e-6);

    % Summarize results
    Qq = boolnet.Q;
    TP = sum(sum(Qgt & Qq));     % numCorrectlinks
    FP = sum(sum(~Qgt & Qq));    % numWrongLinks
    FN = sum(sum(Qgt & ~Qq));    % numMissedLinks
    TN = sum(sum(~Qgt & ~Qq));
    TPR = TP / (TP + FN);        % 1 - Prec = numWrongLinks / numTotalLinksRes
    Prec = TP / (TP + FP);       % 1 - TPR =  numMissedLinks / numTotalLinksGT

    netResults(1, indexDataID) = TP;
    netResults(2, indexDataID) = FP;
    netResults(3, indexDataID) = FN;
    netResults(4, indexDataID) = TN;
    netResults(5, indexDataID) = solverOptions.lambda;

    fprintf('DataSet %02d: Type I  Error = %f%% \n', ...
            dataIDList(indexDataID),  (1-Prec)*100);
    fprintf('            Type II Error = %f%% \n\n', (1-TPR)*100);

    %% Write Results in TSV
    fid = fopen(txtName, 'a');
    fspec = '%d\t%d\t%d\t%d\t%.2e\n';
    fprintf(fid, fspec, netResults(:,indexDataID));
    fclose(fid);

    %% Data Visualization
    fig_hl = figure('visible', 'off');
    pg = plot(netobj);
    xlabeltext = sprintf(['Ts = ' num2str(Ts) '(second); ' ...
                        'SNR = ' num2str(SNR) 'dB;\n'...
                        'Threshold of Zeros: ' num2str(netobj.ThresholdZero) '.']);
    xlabel(xlabeltext);
    solverName = SolverNameMap(solverOptions.solver);
    title([solverName.fullName ' (\lambda = ' num2str(solverOptions.lambda) ')']);
    BoolNetWrong = ~Qgt & Qq;
    [s, t] = find(BoolNetWrong');
    highlight(pg, s', t', 'Edgecolor', 'r');

    % Export Figures in PDF
    figName = [resPath, fnamePrefix, '.pdf'];
    set(fig_hl,'Units','Inches');
    pos = get(fig_hl,'Position');
    set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl, figName, '-dpdf', '-r0')

end


%% Time benchmark and Shutdown parpool
netEtime = toc(netTimer);
fprintf('Elapsed time (network inference) is %f seconds. \n\n', ...
        netEtime);
% delete(gcp);      % shutdown parpool


%% ARCHIVE: Model Comparison
%
% dsfobj = dynet(netobj.dsf());
% [fitScores, yModelResponse, xInitCond] = dsfobj.compare(data, 'plot');
