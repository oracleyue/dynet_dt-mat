function boolNets = multi_runs(net, data, varargin)
% function boolNets = multi_runs(net, data, ModelOrders, SolverConfigs)
% MULTI_RUNS is to perform network inference under different running
% parameters settings, e.g. model orders, regularization parameters, etc.
%
% Input:
%     Optional parameter Name-Value pairs:
%       "ModelOrders"       3xN cell: each row, {px1 pxp px1} integers OR
%                           3xN integers
%       "Lambdas"           2xN or 1xN double
%
%     Other control options:
%       "Visible"           'on' (default), 'off'
%       "SparsityMethod"    'gil1' (default), 'sbl'
%       "SysIdMethod"    'arx' (default), 'narx'
%
% Output:
%     boolNets: 1xN cell, each is a structure of bool network, with members:
%       boolNets{i}.Q, boolNets{i}.P
%
% Copyright [2017] <oracleyue>



% ------ BEGIN OF INPUT ARGUMENT PARSING ------

% Optional arguments handling
argParser = inputParser;
argParser.FunctionName = 'multi_runs';

% Default arguments
defaultModelOrders = [1 2 1];
defaultSolverConfig = {'gil1', [0 5e-3], 2};
expectedSysIdMethod = {'arx', 'narx'};


% Adding Name-Value pairs
addRequired(argParser, 'net');
addRequired(argParser, 'data', ...
            @(x) strcmp(class(data), 'iddata'));
addParameter(argParser, 'ModelOrders', defaultModelOrders, @ismatrix);
addParameter(argParser, 'SolverConfigs', defaultSolverConfig, @iscell);
addParameter(argParser, 'SysIdMethod', 'narx', ...
             @(x) any(validatestring(x,expectedSysIdMethod)));

% Parsing input arguments
parse(argParser, net, data, varargin{:});
net = argParser.Results.net;
data = argParser.Results.data;
ModelOrders = argParser.Results.ModelOrders;
SolverConfigs = argParser.Results.SolverConfigs;
SysIdMethod = argParser.Results.SysIdMethod;

% ------ END OF INPUT ARGUMENT PARSING ------


% Setup config variables
SparsityMethod = SolverConfigs{1};
Lambdas = SolverConfigs{2};
NumSegment = SolverConfigs{3};

% Determine the type of multi_run, i.e. in terms of ModelOrders or Lambdas
if size(ModelOrders, 1) > 1
    runType = 'ModelOrders';
elseif size(Lambdas, 1) > 1
    runType = 'Lambdas';
else
    runType = 'SingleRun';
end

% Branches of different multi_run()
switch runType
  % ---------------
  case 'ModelOrders'

    N = size(ModelOrders, 1);  % number of experimental runs

    % Network resconstruction using different model orders
    boolNets = cell(1,N);
    for i = 1:N
        morders = ModelOrders(i,:);

        if strcmp(SysIdMethod, 'arx')
            net.arx(data, morders, SolverConfigs);
            boolNets{i} = net.BoolNet;
        elseif strcmp(SysIdMethod, 'narx')
            net.narx(data, morders, SolverConfigs);
            boolNets{i} = net.BoolNet;
        else
            error('The chosen SysIdMethod is not supported!')
        end
    end

  % ---------------
  case 'Lambdas'

    N = size(Lambdas, 1);  % number of experimental runs

    % Network resconstruction using different model orders
    boolNets = cell(1,N);
    for i = 1:N
        morders = ModelOrders;
        solver_config = {SparsityMethod Lambdas(i,:) NumSegment};

        if strcmp(SysIdMethod, 'arx')
            net.arx(data, morders, solver_config);
            boolNets{i} = net.BoolNet;
        elseif strcmp(SysIdMethod, 'narx')
            net.narx(data, morders, solver_config);
            boolNets{i} = net.BoolNet;
        else
            error('The chosen SysIdMethod is not supported!')
        end
    end

  % ---------------
  case 'SingleRun'
    warning(['The network reconstruction is performed only once. ' ...
             'Check your "ModelOrders" or "SolverConfigs" settings.'])
    morders = ModelOrders;
    solver_config = SolverConfigs;

    if strcmp(SysIdMethod, 'arx')
        net.arx(data, morders, solver_config);
        boolNets = net.BoolNet;
    elseif strcmp(SysIdMethod, 'narx')
        net.narx(data, morders, solver_config);
        boolNets = net.BoolNet;
    else
        error('The chosen SysIdMethod is not supported!')
    end

end
