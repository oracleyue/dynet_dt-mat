% A wrapper to run scripts with multiple settings and datasets in batch.

% Copyright [2017] <oracleyue>
% Last modified on 23 Aug 2017


close all; clear all
run('../../dynet_startup.m')

% select task
task = 'multiple SNRs';
% task = 'multiple #nodes'

% solvers
solverList = {'gil1', 'gsbl'};
solverSelected = 1:2;

% datasets
dataIDList = 1:50;

% setup environments for different tasks
switch task
  case 'multiple SNRs'
    % path and file settings
    dataPath = '../../SimData/MCDataSNRs/dataMAT/';
    resPath  = '../../Results/resSNRs/';
    figName = 'snr_err';

    % runtime settings
    lambdaList = [1.2e-3, 5.5e-3 0
                  1e-4, 9e-3, 1.8e-1];
    SNRList = [80 60 40 20];

    for indexSolver = 1:length(solverSelected)
        % solver
        solver = solverList{solverSelected(indexSolver)};

        % common settings
        p = 40;      % dim of y
        m = 1;       % dim of u

        for index = 1:length(SNRList)
            % solver settings
            lambda = lambdaList(indexSolver, index);
            solverOptions = solverInit(solver, lambda, 0);

            % dataset settings
            SNR = SNRList(index);

            % run the script
            run('./sDynetHeteroMC.m');
        end
    end

  % ------------------------------------------------------

  case 'multiple #nodes'
    % path and file settings
    dataPath = '../../SimData/MCDataNumNodes/dataMAT/';
    resPath  = '../../Results/resNodes/';
    figName = 'node_err';

    % runtime settings
    lambdaList = [1e-4, 9e-3, 1.8e-1;
                  1.2e-3, 5.5e-3 0];
    NodeList = [80 40 20 10];

    for indexSolver = 1:length(solverSelected)
        % solver
        solver = solverList{solverSelected(indexSolver)};

        % common settings
        m   = 1;     % dim of u
        SNR = 60;    % SNR

        for index = 1:length(NodeList)
            % solver settings
            lambda = lambdaList(indexSolver, index);
            solverOptions = solverInit(solver, lambda, 0);

            % dataset settings
            p = NodeList(index);

            % run the script
            run('./sDynetHeteroMC.m');
        end
    end

end % END: switch task

%% Figure plot and export
run('sStatFigurePlot.m')
