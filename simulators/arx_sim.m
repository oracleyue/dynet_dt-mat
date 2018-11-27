function arx_sim(dataOpt, netOpt, signalOpt, show_plot)
% ARX_SIM uses arx_gen.m to generates random ARX models and performs
% simulations.
%
% Inputs:
%   - dataOpt: struct
%       - numData: integer; length of time series
%       - numExpr: integer; number of experiments, i.e. #datasets
%       - idxStart: all points before this time instant discarded
%       - fileName: string; the file name to save data
%   - netOpt: struct
%       - numNodes: number of nodes in dynet; dim of y
%       - sparsity: float in (0,1); sparsity density of dynet
%       - Pb: Boolean matrix that shows #inputs and which nodes to perturb
%       - orderMax: 1x2 integer vector; max order of numerators and
%       denominators;
%       - Ts: float; sampling period
%   - signalOpt: struct
%       - inputType: 'rgs' or 'step'
%       - inputVariance: positive float of Gaussian input signals
%       - SNR: float in (0.1); signal-noise-ratio
%   - show_plot: Bool; plot example IO signals
% Outputs:
%   export the input/output in the .mat file named "ARXSimData.mat".

% Copyright [2018] <oracleyue>
% Last modified on 19 Jul 2018


% paths
% addpath('./members')    % supportive functions

% parse arguments
if nargin < 4
    show_plot = 0;
end

if nargin < 1
    % data
    num_ts = 1000;          % #points of time series
    num_exp = 1;            % #experiments or #models
    idx_start = 300;        % discard points in signals before idx_start
    fileName = 'ARXSimData.mat';  % file name to save data
    % network
    num_nodes = 10;         % #nodes
    spar = .2;              % density of network sparsity
    Pb = [1; zeros(num_nodes-1,1)];  % boolean P matrix
    den_order_max = 2;      % max order of denominators
    num_order_max = 1;      % max order of numerators
    order_max = [num_order_max den_order_max];
    Ts = .1;                % sampling period
    % signals
    inputType = 'gauss';
    inputVar = 1;           % input signals: multivariate Gaussian
    snr = 1;                % unit: dB
end

if nargin >= 1 && nargin < 3
    error(['You need to specify "dataOpt, netOpt, signalOpt" or use ' ...
           'empty. Read the help doc. '])
end

if nargin >= 3
    % data
    num_ts = dataOpt.numData;
    num_exp = dataOpt.numExpr;
    idx_start = dataOpt.idxStart;
    fileName = dataOpt.fileName;
    % network
    num_nodes = netOpt.numNodes;
    spar = netOpt.sparity;
    Pb = netOpt.Pb;
    num_inputs = size(Pb,2);
    order_max = netOpt.orderMax;
    Ts = netOpt.Ts;
    % signals
    inputType = signalOpt.inputType;
    inputVar = signalOpt.inputVariance;
    snr = signalOpt.SNR;
end

% reparse certain variables
num_ts = num_ts + idx_start;
num_inputs = size(Pb,2);
noiseVar = inputVar / 10^(snr/10);

% storage variables
inputs = cell(num_exp,1);
outputs = cell(num_exp,1);
% models  : struct; save models as ground truth
% simInfo : struct; save simulation configs

% generate models
for idx = 1:num_exp
    [A, B, Qb] = arx_gen(num_nodes, spar, order_max, Pb);
    for k = 1:num_nodes
        C{k,1} = 1;
    end
    sys = idpoly(A,B,C,[],[], noiseVar*eye(num_nodes), Ts);

    switch inputType
      case 'gauss'
        u = idinput([num_ts size(Pb,2)], 'rgs', [0,inputVar]);
      case 'step'
        u = ones(num_ts, size(Pb,2))*inputVar;
    end
    % noise = randn(num_ts, size(Qb,1)) * sqrt(noiseVar);
    opt = simOptions('AddNoise', true);
    y = sim(sys, u, opt);
    inputs{idx} = u(idx_start+1:end,:);
    outputs{idx} = y(idx_start+1:end,:);

    models.A{idx} = A;
    models.B{idx} = B;
    models.Qb{idx} = Qb;
end

% show examples of IO signals
if show_plot
    % idxExp = randi(num_exp);
    for idxExp = 1:num_exp
        if strcmp(fileName, 'ARXSimData.mat')
            figure;
        else
            fighl = figure('visible', 'off');
        end
        % plot time series
        time = 1:size(inputs{idxExp},1);
        subplot(4,1,1)
        idxChannel = randi(num_inputs);
        plot(time, inputs{idxExp}(:,idxChannel), '-');
        ylabel(sprintf('input(%u)',idxChannel))
        title('Examples of input and output signals')
        subplot(4,1,2)
        plot(time, outputs{idxExp}(:,1), '-');
        ylabel('output(1)')
        subplot(4,1,3)
        idxChannel = randi(num_nodes);
        plot(time, outputs{idxExp}(:,idxChannel), '-');
        ylabel(['output(' num2str(idxChannel) ')'])
        subplot(4,1,4)
        plot(time, outputs{idxExp}(:,end), '-');
        ylabel(sprintf('output(%u)', num_nodes))
        xlabel(sprintf('time (Ts = %g)', Ts))
        % export into pdf
        if ~strcmp(fileName, 'ARXSimData.mat')
            figname = [strrep(fileName, '.mat', '_No') ...
                       num2str(idxExp), '.pdf'];
            set(fighl,'Units','Inches');
            pos = get(fighl,'Position');
            set(fighl,'PaperPositionMode','Auto', ...
                      'PaperSize',[pos(3), pos(4)]);
            print(fighl, figname, '-dpdf', '-r0')
        end
    end
end

% export data in .mat
simInfo.SNR = snr;
simInfo.Ts = Ts;
simInfo.inputVariance = inputVar;
save(fileName, 'inputs', 'outputs', 'models', 'simInfo');