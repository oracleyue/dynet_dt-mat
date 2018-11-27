function [cOutputs, input, Ts, boolnet] = ...
    dynet_sim_hetero(model, varargin)
% DYNET_SIM simulates the continuous-time system to generate data for system
% identification.
%
% INPUTS:
%   model: structure
%      - n: integer
%           the dimension of state variables (A-matrix)
%      - p: integer
%           the dimension of output variables
%      - m: integer
%           the dimension of input variables
%      - N: integer
%           the length of time series
%      - IdxNodes: integer vector of dim "p"
%           the indices of stable variables that are outputs;
%           to construct C-matrix
%      - TsRatio: double (recommend > 40)
%           the sampling period = Ts (system aliasing) / TsRatio
%      - sigma: double
%           standard deviation of nosie, i.e. the square-root of variance
%      - x0: scalar or vector
%           scalar: initial states = x0 * randn(n,1)
%           vector: initial states = x0
%      - SparseDensity : (0,1]
%           the sparsity density of A-matrix
%   numdata: integer; (default: 1, i.e. no heterogeneous datasets)
%   optional parameter Name-Value pairs:
%      "PathName": string; folder to save data and figures
%      "FileName": string; names of data and figures
%      "Input": structure
%         - type: string; "step", "gauss"
%         - param: mx1 or 1x1 vector ("step"); 1x2 matrix ("gauss")
%                  "step": magnitude
%                  "gauss": mean and variance;
%      "Mode": string; control the way to construct B,K
%         - 'simple': input and noise are fed in at nodes
%         - 'random': randomly choose states to fed in input and noise
%      "Export": {-1, 0, 1, 2, 3}
%         -1 : inhibits the exportation of data and figures
%          0 (default): save data and figures; NOT show figures
%          1 : save data and plot/save/show figures
%          2 : save data; NOT plot/save/show figures
%          3 : plot/save/show figures; NOT save data
%          4 : plot/show figures; but NOT save figures or data
%
% OUTPUTS:
%   input: Nxm matrix
%   output:
%     - iddata structure, if using "sysid_toolbox" for simulation;
%       ("output.y" for Nxp output signals, "output.Ts" for Ts)
%     - structure, if using "euler-maruyama" for simulation;
%       ("output.y" for Nxp output signals, "output.Ts" for Ts)
%   boolnet: structure
%     members: "Q", "P"; the ground truth boolean network
%
% FILE/EXPORT OPERATIONS (unless "Export" = -1):
%     - create folders "dataCSV", "dataMAT", "figures"
%     - save data and save/plot figures (see "Export")

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 23 Aug 2017



% --- BEGIN OF INPUT ARGUMENT PARSING ---
argParser = inputParser;
argParser.FunctionName = 'dynet_sim_hetero';

% default optional parameters
defaultNumDatasets = 1;
defaultPathName = './SimData/';
defaultFileName = 'dynetsim_test';
defaultInput.type = 'step';
defaultInput.param = 1;
defaultMode = 'simple';
defaultExport = 0;
validExport = [-1, 0, 1, 2, 3, 4];

% adding name-value pairs
addRequired(argParser, 'model', @isstruct);
addOptional(argParser, 'numdata', defaultNumDatasets, ...
            @(x) isscalar(x) && rem(x,1) == 0);
addParameter(argParser, 'PathName', defaultPathName, @ischar);
addParameter(argParser, 'FileName', defaultFileName, @ischar);
addParameter(argParser, 'Input', defaultInput, @isstruct);
addParameter(argParser, 'Mode', defaultMode, ...
             @(x) ischar(x) && (strcmp(x, 'simple') || ...
                                strcmp(x, 'random') || ...
                                strcmp(x, 'measurement') || ...
                                strcmp(x, 'oneprocess')) );
addParameter(argParser, 'Export', defaultExport,...
             @(x) isscalar(x) && any(x == validExport));

% parsing arguments
parse(argParser, model, varargin{:});
model = argParser.Results.model;
numdata = argParser.Results.numdata;
pname = argParser.Results.PathName;
fname = argParser.Results.FileName;
inputStruct = argParser.Results.Input;
mode = argParser.Results.Mode;
exportValue = argParser.Results.Export;
% --- END OF INPUT ARGUMENT PARSING ---


%% Systems Generation and Simulation

% Message on the end of simpulation
fprintf('DataSet ID: "%s"\n', fname);
fprintf('INFO: starting ...\n');

% Unpack basic settings
n = model.n; p = model.p; m = model.m;
N = model.N; Tsr = model.TsRatio;
sigma = model.sigma; InitVal = model.x0;
% IdxNodes = sort(model.IdxNodes);  % sort in case of random orders
IdxNodes = model.IdxNodes;
SparseDensity = model.SparseDensity;

% Setup state-space models for simulation

% generating sparse stable A
% [A, fig_hl_A] = dynet_sprandstab(n, 'nicolo', .8, [5 4]);
[A, fig_hl_A] = dynet_sprandstab(n, 'nicolo-overlap', SparseDensity, 4);

% determine the sampling frequency (= SystemAliasingFrequency / TsRatio)
eigA = eig(A);
w_sysali = 2*max(abs(imag(eigA)));
Ts_sysali = 2*pi/w_sysali;
Ts = Ts_sysali/Tsr;

% additional arguments of state-space models
C = eye(n);
idxC = IdxNodes;
C = C(idxC,:);
D = zeros(p, m);
switch mode
  case 'simple'
    if m > p
        error(['The "simple" mode only supports "dim input <= dim output". Change ' ...
               'to the "random" mode.']);
    end
    B = eye(n);
    idxB = idxC(1:m);
    B = B(:,idxB);   % input is fed in at output nodes, instead of hidden states

    K = eye(n);
    K = K(:, idxC); % process noises added to output nodes

  case 'random'
    B = eye(n);
    randIndexB = randperm(n);
    idxB = sort(randIndexB(1:m)); % random states to feed inputs
    B = B(:,idxB);

    K = eye(n);
%     randIndexK = randperm(n);
%     idxK = sort(randIndexK(1:p));
%     K = K(:, idxK); % process noises added to random p nodes
    K = K(:, idxC); % process noises added to output nodes    

  case 'measurement'  % FOR DEBUG ONLY: no process noise
    B = eye(n);
    idxB = idxC(1:m);
    B = B(:,idxB);   % input is fed in at output nodes

    K = zeros(n);
    K = K(:, idxC); % process noises added to output nodes

  case 'oneprocess'  % FOR DEBUG ONLY: one channel of process noise +
                     % measurement noise
    B = eye(n);
    idxB = idxC(1:m);
    B = B(:,idxB);   % input is fed in at output nodes

    col = randi(length(idxC),1);  % randomly pick one output channel
    Kcol = eye(n);
    Kcol = Kcol(:,idxC(col));
    K = zeros(n, length(idxC));
    K(:,idxC(col)) = Kcol; % process noise added to one output node
end

% initial values
if isscalar(InitVal)
    init_val = randn(n,1)*InitVal;
else
    if size(InitVal) == [n,1]
        init_val = InitVal;
    else
        error('model.x0 should be either a scalar or nx1 vector.')
    end
end
% input signals
switch inputStruct.type
  case 'step'
    if isscalar(inputStruct.param)
        input = [zeros(round(N/10), m); ones(N-round(N/10), m)]...
                * inputStruct.param;
    else
        input = bsxfun(@times, ...
                       [zeros(round(N/10), m); ones(N-round(N/10), m)], ...
                       inputStruct.param');
    end
  case 'gauss'
    input = idinput([N m], 'rgs', [0,1], ...
                    [inputStruct.param(1) - inputStruct.param(2), ...
                     inputStruct.param(1) + inputStruct.param(2)]);
end
% noise signals
noise = idinput([N,p], 'rgs', [0,1], [-sigma, sigma]);  % use sysid toolbox
% noise = randn(N, p)*sigma;  % if N is small, mean is not very close to 0

% Simulating State-space models
cOutputs = {};
cAs = {};
for k = 1:numdata
    % generate heterogeneous systems
    if k == 1,  Asim = A;
    else
        Asim = pertubeAmatrix(A);
        while ~isStable(Asim)
            Asim = pertubeAmatrix(A);
        end
    end
    cAs = [cAs Asim];
    fprintf('      ...sparse stable A generated\n');

    % perform simulation
    method = 'sysid_toolbox';
    switch method
      case 'sysid_toolbox'
        % --- simulations using sysid toolbox ---
        sys_ss = idss(Asim, B, C, D, K, 'Ts', 0);
        simopt = simOptions('AddNoise', true, 'NoiseData', noise, ...
                            'InitialCondition', init_val);
        output = sim(sys_ss, iddata([], input, Ts), simopt);

        cOutputs = [cOutputs output.y]; % save signals in cells

      case 'euler-maruyama'
        % ---  simulations using Euler-Maruyama algorithm ---
        sys.A = Asim; sys.B = B; sys.C = C; sys.D = D; sys.K = K;
        data.u = input; data.Ts = Ts; data.sigma = sigma;
        data.X0 = init_val; data.N = N;

        output.y = ss_sim(sys, data);
        output.Ts = Ts;

        cOutputs = [cOutputs output.y]; % save signals in cells
    end
end


%% Dynamic Network Models: calculate QP from A,B

% change orders of nodes in A such that C = [I 0], i.e. x= [y z]'
idxPerm = [IdxNodes setdiff(1:1:n, IdxNodes)];
P = eye(n); P = P(idxPerm,:);
A = P*A*inv(P); B = P*B;
% calculate DSF
[Q, P] =getpq(A, B, p);
% find boolean network
s = 1;    Q1 = double(subs(Q)); P1 = double(subs(P));
s = .11;  Q2 = double(subs(Q)); P2 = double(subs(P));
s = .111; Q3 = double(subs(Q)); P3 = double(subs(P));
Qgt = logical(Q1) | logical(Q2) | logical(Q3);
Pgt = logical(P1) | logical(P2) | logical(P3);
% structure for output arguments
DSF.Q = Q;
DSF.P = P;
boolnet.Q = Qgt;
boolnet.P = Pgt;


%% Exporting data and figures

% Flag Interpretations
exportFlag = ~(exportValue == -1);
dataFlag = (exportValue == 0 || exportValue == 1 || exportValue == 2);
plotFlag = ~(exportValue == 2);
showFlag = plotFlag & ~(exportValue == 0);
pdfFlag = plotFlag & (exportValue == 0 || exportValue == 1 || exportValue == 3);

% ---- BEGIN: EXPORTING ---
if exportFlag

% Plot Figures
if plotFlag
    % % plot output signals
    % fig_hl = figure('visible', 'off');
    % % set(fig_hl,'Units','Inches', 'position',[2.1528 3.2639 12 21]);
    % set(fig_hl,'Units','Inches');
    % time = (1:N)*Ts;
    % % plot output signals in subplot
    % for i = 1:p
    %     subplot(ceil(p/4),4,i)
    %     plot(time, output.y(:,i), '-')
    %     xlabel('time')
    %     ylabel(['y_{' num2str(i) '}'])
    % end

    % plot output signals in one figure for preview
    fig_hl_all = figure('visible', 'off');
    set(fig_hl_all,'Units','Inches');
    plot((1:N)', cOutputs{1}, '-')
    xlabel('index'); ylabel('outputs')

    % plot diagraph
    fig_hl_Q = figure('visible', 'off');
    set(fig_hl_Q,'Units','Inches', 'position',[9.4306 7.5000 7.7778 6.0694]);
    % [target, source] = find(Qgt);
    % Ggt = digraph(source, target);
    Ggt = digraph(Qgt);
    pg = plot(Ggt, 'Layout', 'force');
    pg.NodeColor = 'red';
    title('Ground truth of Boolean network')
end
% Show Figures
if showFlag
    set(fig_hl, 'visible', 'on');
    set(fig_hl_all, 'visible', 'on');
    set(fig_hl_A, 'visible', 'on');
    set(fig_hl_Q, 'visible', 'on');
end


% Data Export
if dataFlag

    % variable for sim info
    ss_dim = n;      % dimension of ss sys, i.e. A
    out_dim = p;
    in_dim = m;
    ts_dim = N;      % length of time series
    Ts_sim = Ts;
    x0_sim = init_val;

    % %
    % % Export Datasets in CSV Format
    % %
    % if ~exist([pname 'dataCSV'], 'dir')
    %     mkdir([pname 'dataCSV']);
    % end

    % % saving sim info in CSV files
    % txtfname_info = [pname 'dataCSV/' fname '_info.txt'];
    % fid = fopen(txtfname_info, 'w');
    % % Qgt:
    % Format1='%f ';
    % for i=1:p-1, Format1=[Format1,'\t %f']; end
    % Format1=[Format1,'\n'];
    % fprintf(fid,Format1, Qgt');
    % % Ts:
    % fprintf(fid,'\n');
    % Format2='%f\n';
    % fprintf(fid,Format2, Ts);
    % % A:
    % fprintf(fid,'\n');
    % Format3='%f ';
    % for i=1:n-1, Format3=[Format3,'\t %f']; end
    % Format3=[Format3,'\n'];
    % fprintf(fid,Format3, A');
    % % w_sysali:
    % fprintf(fid,'\n');
    % Format4='%f\n';
    % fprintf(fid,Format4, w_sysali);
    % fclose(fid);

    % % saving sim signals in CSV files
    % txtfname_out = [pname 'dataCSV/' fname '_out.txt'];
    % txtfname_in = [pname 'dataCSV/' fname '_in.txt'];
    % % saving output signals
    % fid = fopen(txtfname_out,'w');
    % Format='%f ';
    % for i=1:p-1, Format=[Format,'\t %f']; end
    % Format=[Format,'\n'];
    % fprintf(fid,Format,output.y');
    % fclose(fid);
    % % saving input signals
    % fid = fopen(txtfname_in,'w');
    % Format='%f ';
    % for i=1:m-1, Format=[Format,'\t %f']; end
    % Format=[Format,'\n'];
    % fprintf(fid,Format,input');
    % fclose(fid);

    %
    % Export Datasets in ".mat" Format
    %
    if ~exist([pname 'dataMAT'], 'dir')
        mkdir([pname 'dataMAT']);
    end
    % saving in mat format
    %   note: save/load not available inside /parfor/
    matfname = [pname 'dataMAT/' fname '_info.mat'];
    save(matfname, 'ss_dim', 'out_dim', 'in_dim', 'ts_dim', ...
         'cAs', 'B', 'C', 'D', 'K', ...
         'Qgt', 'IdxNodes', ...
         'sigma', 'Ts_sim', 'x0_sim', 'noise', 'w_sysali');
    matfname = [pname 'dataMAT/' fname '.mat'];
    save(matfname, 'Ts', 'input', 'cOutputs');
end

% Save Figures in PDF
if pdfFlag
    % create the folder "figures" to save figs
    if ~exist([pname 'figures'], 'dir')
        mkdir([pname 'figures']);
    end

    % figname = [pname 'figures/' fname '_output.pdf'];
    % pos = get(fig_hl,'Position');
    % set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print(fig_hl, figname,'-dpdf','-r0')

    figname = [pname 'figures/' fname '_output_preview.pdf'];
    pos = get(fig_hl_all,'Position');
    set(fig_hl_all,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl_all, figname,'-dpdf','-r0')

    figname = [pname 'figures/' fname '_A.pdf'];
    pos = get(fig_hl_A,'Position');
    set(fig_hl_A,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl_A, figname,'-dpdf','-r0')

    figname = [pname 'figures/' fname '_Q.pdf'];
    pos = get(fig_hl_Q,'Position');
    set(fig_hl_Q,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_hl_Q, figname,'-dpdf','-r0')
end

end
% ---- END: EXPORTING ---

% Message on the end of simpulation
fprintf('      simulation finished.\n\n');

end




% ================================================================
% Local Functions
% ================================================================

%
% check stability of A-matrix
%
function res = isStable(A)
eig_A = eig(A);
res = ~sum(real(eig_A) >= 0);
end

%
% perturb a few parameters (off-diagonal) in A-matrix
%
function A = pertubeAmatrix(A)
A_nodiag = A - diag(diag(A));  % remove diagonal elements

vIndex = find(A_nodiag ~= 0);
minVal = min(min(abs(A_nodiag)));
perbVector = sprandn(length(vIndex),1, .3)*minVal*.01;
perbVector = double(perbVector);

A(vIndex) = A(vIndex) + perbVector;
end
