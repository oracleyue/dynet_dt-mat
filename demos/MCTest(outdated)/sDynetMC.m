% Perform network inference on multiple datasets.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 04 Aug 2017



close all
clear all

%% Environment Settings

% path settings
run('../../dynet_startup.m')
dataPath = '../../SimData/MCDataSNRs/dataCSV/';
resPath = '../../Results/';

% solver settings
global glb_optimizer
glb_optimizer = 'cvx';
global MAX_ITER ABSTOL  %precision for prox
MAX_ITER = 1e4;
ABSTOL = 1e-8;
solver = 'gil1';  lambda = 1e-4;
% solver = 'gsbl';  lambda = 2e-3*4;

% dataset settings
p = 40;             % dim of y
SNR = 60;           % signal-to-noise ratio
m = 1;              % dim of u
dataIDList = 1:1;   % list of dataID

% variables: memory, I/O, files
netResults = zeros(5, length(dataIDList));
solverName = SolverNameMap(solver);
txtName = [resPath 'dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_' ...
           solverName.abbrName '.txt'];
% if exist(txtName, 'file'), delete(txtName); end
fid = fopen(txtName, 'w');
fspec = 'TP\tFP\tFN\tTN\tlambda\n';  % print the column names for tables
fprintf(fid, fspec); fclose(fid);



%% Network Inference

% parpool('local', 4);
netTimer = tic;

for indexDataID = 1:length(dataIDList)
    dataID = dataIDList(indexDataID);

    % Load Datasets
    fnamePrefix = ['dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_Data'...
                   sprintf('%02d',dataID)];
    fname = [dataPath fnamePrefix '_out.txt'];
    output = dlmread(fname, '\t');
    fname = [dataPath fnamePrefix '_in.txt'];
    input = dlmread(fname, '\t');
    fname = [dataPath fnamePrefix '_info.txt'];
    simInfo = dlmread(fname, '\t');
    Qgt = simInfo(1:p, 1:p);
    Ts_sim = simInfo(p+1,1);

    % Build iddata
    Nsim = size(output, 1);  % time series length
    DT = 1;   % sampling freq: Ts = Ts_sim * DT
    Ts = Ts_sim * DT;
    output = output(100:DT:Nsim, :);
    input = input(100:DT:Nsim, :);
    data = iddata(output, input, Ts);

    % Call network inference function
    netobj = idnet(p,m);
    norders = [2 5 5];
    netobj.ThresholdZero = 1e-6; % used for generating bool network
    netobj.arx(data, norders, {solver, lambda});
    Qq = logical(netobj.BoolNet.Q);

    % Summarize inference results
    numTotalLinksGT = sum(sum(Qgt));
    numTotalLinksRes = sum(sum(Qq));
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
    netResults(5, indexDataID) = lambda;

    fprintf('Data Set %02d: Type I  Error = %f%% \n', ...
            dataIDList(indexDataID),  (1-Prec)*100);
    fprintf('             Type II Error = %f%% \n\n', (1-TPR)*100);


    % Write Results in TSV
    fid = fopen(txtName, 'a');
    fspec = '%f\t%f\t%f\t%f\t%f\n';
    fprintf(fid, fspec, netResults(:,indexDataID)');
    fclose(fid);


    %% Data Visualization and Export


    % Plot digraphs
    fig_hl = figure('visible', 'off');
    pg = plot(netobj);
    xlabeltext = sprintf(['Sampling period = ' num2str(DT) 'Ts (second);\n ' ...
                        'Samples: 1:' num2str(DT) ':' num2str(Nsim) '; '...
                        'Signal-to-Noise Ratio: ' num2str(SNR) ';\n'...
                        'Threshold of zero: ' num2str(netobj.ThresholdZero) '.']);
    xlabel(xlabeltext)
    title([solverName.fullName ' (\lambda = ' num2str(lambda) ')']);
    % highlight wrong links in red
    BoolNetWrong = ~Qgt & Qq;
    [s, t] = find(BoolNetWrong');
    highlight(pg, s', t', 'Edgecolor', 'r');

    % Export figures in PDF
    figName = [resPath, fnamePrefix, '_', solverName.abbrName, '.pdf'];
    set(fig_hl,'Units','Inches');
    pos = get(fig_hl,'Position');
    set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl, figName, '-dpdf', '-r0')

end   % END (indexDataID): network inference of multiple datasets

% time benchmark and close
netEtime = toc(netTimer);
fprintf('The elapsed time of network inference: %f seconds \n', ...
        netEtime);
% delete(gcp);      % shutdown parpool





%% ARCHIVE: Model Comparison
%
% dsfobj = dynet(netobj.dsf());
% [fitScores, yModelResponse, xInitCond] = dsfobj.compare(data, 'plot');
