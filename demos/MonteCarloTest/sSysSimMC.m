% This is a script to simulate continuous-time DSF models to generate
% all datasets for the Automatica paper.

% INFO: It simulates the state-space systems with random sparse stable
% A-matrices.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 03 Aug 2017 14:41:18



% Global Workspace Settings
% --------------------------------------------------------------
% gSelectTask = 'multiple #nodes';
gSelectTask = 'multiple SNRs';
gNumDataSets = 2;
gNodeList = [10 20 40 80];
gSNRList = [60 40 20 0]; %SNR = 10*log10(var_u / var_e)
gRGSsigma = 1; % sigma when Gaussian; magnitude when step

% Log Files: create or clean
switch gSelectTask
  case 'multiple #nodes'
    pathName = '../../SimData/MCDataNumNodes/';
  case 'multiple SNRs'
    pathName = '../../SimData/MCDataSNRs/';
end
fid = fopen([pathName 'dataList.log'], 'w');
fprintf(fid, 'List of data sets:\n------\n');
fclose(fid);

% Start timer
simTimer = tic;

% Generate Datasets in Batch
for indexJob = 1:gNumDataSets
    switch gSelectTask

        % Experiment 1: multiple #nodes
        % usage: To perform it, set /gSelectTask/
        % --------------------------------------------------------------
      case 'multiple #nodes'
        for indexNode = 1 %1:length(gNodeList)
            % Model Setup
            SNR = 60;
            p = gNodeList(indexNode);

            model.p = p;                            % dim of outputs
            model.n = model.p * 3;                  % dim of states

            % - choose hidden states: fixed equal number in each edge
        %     numHiddenNodes = 1;
        %     model.IdxNodes = 1:(numHiddenNodes+2):model.n;
            % - choose hidden states: randomly select
        %     indexPerm = randperm(model.n);
        %     model.IdxNodes = indexPerm(1:model.p);
            % - choose hidden states: randomly select; >= 1 in each edge
            indexPerm = randperm(model.n);
            indexPermOddIndices = find(rem(indexPerm,2) ~= 0);
            indexPermOdd = indexPerm(indexPermOddIndices);
            model.IdxNodes = indexPermOdd(1:model.p);

            model.m = 1;                      % dim of inputs
            model.N = 1000;                         % #samples
            model.TsRatio = 60;                     % see doc of "dynet_sim.m"
            model.sigma = sqrt(gRGSsigma^2 / 10^(SNR/10));
            model.x0 = 0;                           % initial value of states
            model.SparseDensity = .002/(gNodeList(indexNode)/10);

            % System Simulation
            inputSignal.type = 'gauss';
            inputSignal.param = [0, gRGSsigma];
            %     inputSignal.type = 'step';
            %     inputSignal.param = gRGSsigma;
            pathName = '../../SimData/MCDataNumNodes/'; %  path to save data/figures
            fileName = ['dynet' '_p' num2str(model.p) '_SNR' num2str(SNR) 'dB'...
                        '_Data' sprintf('%02d', indexJob) ];
            dynet_sim(model, 'PathName', pathName, 'FileName', fileName, ...
                      'Input', inputSignal, 'Mode', 'measurement', ...
                      'Export', 0);

            % Log of datasets
            fid = fopen([pathName 'dataList.log'], 'a');
            fprintf(fid, [fileName '\n']);
            fclose(fid);
        end


        % Experiment 2: multiple SNRs
        % usage: To perform it, set /gSelectTask/
        % --------------------------------------------------------------
      case 'multiple SNRs'
        for indexSNR = 1 %:length(gSNRList)
            % Model Setup
            p = 40;
            SNR = gSNRList(indexSNR);

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

            model.m = 1;                            % dim of inputs
            model.N = 1000;                         % #samples
            model.TsRatio = 40;                     % see doc of "dynet_sim.m"
            model.sigma = sqrt(gRGSsigma^2 / 10^(SNR/10));
            model.x0 = 0;                           % initial value of states
            model.SparseDensity = .02/p;

            % System Simulation
            inputSignal.type = 'gauss';
            inputSignal.param = [0, gRGSsigma];
%             inputSignal.type = 'step';
%             inputSignal.param = gRGSsigma;
            pathName = '../../SimData/MCDataSNRs/'; %  path to save data/figures
            fileName = ['dynet' '_p' num2str(model.p) '_SNR' num2str(SNR) 'dB'...
                        '_Data' sprintf('%02d', indexJob) ];
            % dynet_sim(model, 'PathName', pathName, 'FileName', fileName, ...
            %           'Input', inputSignal, 'Mode', 'random', ...
            %           'Export', 0);
            dynet_sim_hetero(model, 2, 'PathName', pathName, 'FileName', fileName, ...
                             'Input', inputSignal, 'Mode', 'random', ...
                             'Export', 0);

            % Log of datasets
            fid = fopen([pathName 'dataList.log'], 'a');
            fprintf(fid, [fileName '\n']);
            fclose(fid);
        end

    end % END-SWITCH: gSelectTask
end % END-FOR: indexjob

% End timer
simEtime = toc(simTimer);
fprintf('Elapsed time (model simulation) is %f seconds \n\n', ...
        simEtime);
