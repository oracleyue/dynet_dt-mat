% This script uses "arx_sim.m" (call "arx_gen" inside) to randomly
% generate sparse stable ARX models.
% Info:
%   - run "dynet_startup.m" first to setup search paths;
%   - "arx_sim.m": ProjectRoot/simulaters/.

% Copyright [2018] <oracleyue>
% Last modified on 19 Jul 2018


type = 3;

switch type
  case 1
    % use default configs and save signals in "ARXSimData.mat"
    arx_sim();

  case 2
    % quick generate configs and save signals in "mydata.mat"
    [dataOpt, netOpt, signalOpt] = arxsimOptions();
    dataOpt.fileName = 'mydata.mat';
    show_plot = 1;   % optional; plot example signals
    arx_sim(dataOpt, netOpt, signalOpt, show_plot);

  case 3
    % modify more configs (refer to help doc of "arx_sim.m")
    [dataOpt, netOpt, signalOpt] = arxsimOptions();
    dataOpt.numData = 1000;
    dataOpt.numExpr = 5;        % 5 random ARXs and hence 5 datasets
    netOpt.numNodes = 8;        % dim of network
    netOpt.sparsity = .3;       % sparsity density
    signalOpt.SNR = 10;         % SNR in dB
    arx_sim(dataOpt, netOpt, signalOpt);
end
