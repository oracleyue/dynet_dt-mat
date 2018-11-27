function [dataOpt, netOpt, signalOpt] = arxsimOptions(p,m)
% ARXSIMOPTIONS initializes the structures used for arx_sim.m, for
% quick usages.

% Copyright [2018] <oracleyue>
% Last modified on 19 Jul 2018


% parse arguments
if nargin == 0
    p = 10; m = 1;
end

% data
dataOpt.numData = 1000;        % #points of time series
dataOpt.numExpr = 1;           % #experiments or #models
dataOpt.idxStart = 0;          % discard points before this time
dataOpt.fileName = 'ARXSimData.mat';     % file name to save data

% network
netOpt.numNodes = p;           % #nodes
netOpt.sparity = .2;           % density of network sparsity
% netOpt.Pb = [1; zeros(netOpt.numNodes-1,1)];
matI = eye(p);
netOpt.Pb = matI(:,m);         % boolean P matrix
netOpt.orderMax = [1 2];       % max orders of num and den polynomials
netOpt.Ts = .1;                % sampling period
% signals
signalOpt.inputType = 'gauss'; % gaussian input; or 'step'
signalOpt.inputVariance = 1;   % inputs: multivariate Gaussian
signalOpt.SNR = 1;             % unit: dB
