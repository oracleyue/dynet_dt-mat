% Perform network inference on one hetergeneous dataset, which consists
% of input-output time series of multiple experiments.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 23 Aug 2017



close all
clear all

%% Environment Settings

% path settings
run('../../dynet_startup.m')
dataPath = '../../SimData/MCDataSNRs/dataMAT/';
resPath  = '../../Results/';

% solver settings
solverOptions = solverInit('gsbl', 5e-3, 0);
% solverOptions = solverInit('gil1', 8e-4);

% dataset settings
p      = 40;    % dim of y
m      = 1;     % dim of u
SNR    = 60;    % signal-to-noise ratio
dataID = 1;     % dataset ID No

% flags
EXPORT = 1;    % export figures and results


%% Network Inference

% Start parpool and timer
% parpool('local', 20);
netTimer = tic;

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
Prec = TP / (TP + FP);       % 1 - Prec = numWrongLinks / numTotalLinksRes
TPR  = TP / (TP + FN);       % 1 - TPR =  numMissedLinks / numTotalLinksGT
TypeIError = (1-Prec)*100;
TypeIIError = (1-TPR)*100;

fprintf('Type I  Error: %f%% \n', (1-Prec)*100);
fprintf('Type II Error: %f%% \n\n', (1-TPR)*100);


%% Data Visualization and Export

% Plot Networks
fig_hl = figure;  %figure('visible', 'off');
pg = plot(netobj);
xlabeltext = sprintf(['Ts = ' num2str(Ts) '(second); ' ...
                    'SNR = ' num2str(SNR) 'dB;\n'...
                    'Threshold of Zeros: ' num2str(netobj.ThresholdZero) '.']);
xlabel(xlabeltext);
solverName = SolverNameMap(solverOptions.solver);
title([solverName.fullName ' (\lambda = ' num2str(solverOptions.lambda) ')']);
% highlight wrong links in red
BoolNetWrong = ~Qgt & Qq;
[s, t] = find(BoolNetWrong');
highlight(pg, s', t', 'Edgecolor', 'r');

%% Export
if EXPORT
    % Export figures in PDF
    figName = [resPath, fnamePrefix, '.pdf'];
    set(fig_hl,'Units','Inches');
    pos = get(fig_hl,'Position');
    set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl, figName, '-dpdf', '-r0')

    % Write Results in TSV
    netResults(1) = TP;  netResults(2) = FP;
    netResults(3) = FN;  netResults(4) = TN;
    netResults(5) = solverOptions.lambda;
    netResults(6) = TypeIError;
    netResults(7) = TypeIIError;

    txtName = [resPath fnamePrefix '.txt'];    
    if ~exist(txtName, 'file')
        fid = fopen(txtName, 'w');
        fspec = 'TP\tFP\tFN\tTN\tlambda\tTypeI\tTypeII\n';
        fprintf(fid, fspec);
        fclose(fid);
    end
    fid = fopen(txtName, 'a');
    fspec = '%d\t%d\t%d\t%d\t%.2e\t%.2f\t%.2f\n';
    fprintf(fid, fspec, netResults');
    fclose(fid);
end


%% Time benchmark and Shutdown parpool
netEtime = toc(netTimer);
fprintf('Elapsed time (network inference) is %f seconds \n', ...
        netEtime);
% delete(gcp);      % shutdown parpool


%% ARCHIVE: Model Comparison
%
% dsfobj = dynet(netobj.dsf());
% [fitScores, yModelResponse, xInitCond] = dsfobj.compare(data, 'plot');
