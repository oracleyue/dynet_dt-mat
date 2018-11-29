% Perform network inference on a sepecific dataset

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
solverOptions = solverInit('gsbl', 1.5e-3, 1);
% solverOptions = solverInit('gil1', 3e-3);

% dataset settings
p = 40;             % dim of y
m = 1;              % dim of u
SNR = 80;           % signal-to-noise ratio
dataID = 1;         % dataset ID No

% flags
tsvFlag = 0;    % export results in tsv
pdfFalg = 0;    % export figures in pdf


%% Network Inference

% parpool('local', 4);
netTimer = tic;

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
Nsim = size(output, 1);   % time series length
DT = 1;   % sampling freq: Ts = Ts_sim * DT
Ts = Ts_sim * DT;
output = output(1:DT:Nsim, :);
input = input(1:DT:Nsim, :);
data = iddata(output, input, Ts);

% Call network inference function
netobj = idnet(p,m);
norders = [3 3 3];
netobj.ThresholdZero = 1e-6; % used for generating bool network
netobj.arx(data, norders, solverOptions);
% boolnet = netobj.bool('L', 1e-6);
boolnet = netobj.BoolNet;


%% Data Visualization, Summary and Export

% Plot digraphs
fig_hl = figure;
% fig_hl = figure('visible', 'off');
pg = plot(netobj);
xlabeltext = sprintf(['Sampling period = ' num2str(DT) 'Ts (second);\n ' ...
                    'Samples: 1:' num2str(DT) ':' num2str(Nsim) '; '...
                    'Signal-to-Noise Ratio: ' num2str(SNR) ';\n'...
                    'Threshold of zero: ' num2str(netobj.ThresholdZero) '.']);
xlabel(xlabeltext)
solverName = SolverNameMap(solverOptions.solver);
title([solverName.fullName ' (\lambda = ' num2str(solverOptions.lambda) ')']);
% highlight wrong links in red
Qq = logical(netobj.BoolNet.Q);
BoolNetWrong = ~Qgt & Qq;
[s, t] = find(BoolNetWrong');
highlight(pg, s', t', 'Edgecolor', 'r');

% Export figures in PDF
if pdfFalg
    figName = [resPath, fnamePrefix, '.pdf'];
    set(fig_hl,'Units','Inches');
    pos = get(fig_hl,'Position');
    set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl, figName, '-dpdf', '-r0')
end

% Summarize inference results
numTotalLinksGT = sum(sum(Qgt));
numTotalLinksRes = sum(sum(Qq));
TP = sum(sum(Qgt & Qq));     % numCorrectlinks
FP = sum(sum(~Qgt & Qq));    % numWrongLinks
FN = sum(sum(Qgt & ~Qq));    % numMissedLinks
TN = sum(sum(~Qgt & ~Qq));
TPR = TP / (TP + FN);        % 1 - Prec = numWrongLinks / numTotalLinksRes
Prec = TP / (TP + FP);       % 1 - TPR =  numMissedLinks / numTotalLinksGT

netResults(1) = TP;
netResults(2) = FP;
netResults(3) = FN;
netResults(4) = TN;
netResults(5) = solverOptions.lambda;

fprintf('Type I  Error: %f%% \n', (1-Prec)*100);
fprintf('Type II Error: %f%% \n\n', (1-TPR)*100);

% Write Results in TSV
if tsvFlag
    txtName = [resPath fnamePrefix '.txt'];
    fid = fopen(txtName, 'w');
    fspec = '%f \t %f \t %f \t %f \t %f\n';
    fprintf(fid, fspec, netResults');
    fclose(fid);
end


% time benchmark and close
netEtime = toc(netTimer);
fprintf('The elapsed time of network inferce: %f seconds \n', ...
        netEtime);
% delete(gcp);      % shutdown parpool


%% ARCHIVE: Model Comparison
%
% dsfobj = dynet(netobj.dsf());
% [fitScores, yModelResponse, xInitCond] = dsfobj.compare(data, 'plot');
