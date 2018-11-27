% This script is to simulate continuous-time DSF models to generate datasets.
% INFO: It simulates state-space systems with random sparse stable A-matrices.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 04 Aug 2017



% Simluation Settings
% ----------------------------------------------------
% Note: read the doc of "dynet_sim" to see its meanings
model.n = 20;                           % dim of states
model.p = 10;                           % dim of outputs
% numHiddenNodes = 1;                   % equal number of hidden states
% model.IdxNodes = 1:(numHiddenNodes+1):model.n;
indexPerm = randperm(model.n);
model.IdxNodes = indexPerm(1:model.p);  % randomly select hidden states
model.m = 2;                            % dim of inputs
model.N = 1000;                         % #samples
model.TsRatio = 40;     % choose an appropriate Ts up to A; recommend > 40
model.sigma = .01;      % standard deviration of Gaussian noise
model.x0 = .5;          % or nx1 vector that fixes the initial states
model.SparseDensity = .005/4;           % sparsity of A matrix


% System Simulation
% ----------------------------------------------------

% % Use step input signals
% inputSignal.type = 'step';              % the first 10% of N points are zero
% inputSignal.param = 2;                  % same magnitude for all inputs
% inputSignal.param = randn(model.m,1)*2; % different magnitude

% Use random gaussian white noise as input signals
inputSignal.type = 'gauss';
inputSignal.param = [0, .01];             % [mean, std]


% //Demo 1: simulate model, without saving data/figures into files ('export' = -1)
% [input, output, boolnet] = dynet_sim(model, 'Input', inputSignal, 'Export', -1);
% [input, output] = dynet_sim(model, 'Input', inputSignal, 'Export', -1);


% //Demo 2: simulate and saving data/figures
pathName = '../../SimData/';                        %  path to save data/figures
% dynet_sim(model, 'Input', inputSignal);       % use default path './'
dynet_sim(model, 'PathName', pathName, 'FileName', 'test', ...
          'Input', inputSignal, 'Export', 0);