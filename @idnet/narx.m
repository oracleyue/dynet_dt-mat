function narx_etime = narx(net, rawdata, morders, solver_config)
% IDENT.NARX infers the Boolean nonlinear dynamic network
%    via multiple ARX models
%
% Input:
%     rawdata: iddata object
%          - NOTE: rawdata cannot be cx1 cell, which not supports "narx"
%     morders = {na, nby, nbu}: model orders of A_i, B^y_ij, B^u_ij (C = I)
%          - na : px1
%          - nby: pxp, if px1,then all p columns choosen to be same
%          - nbu: pxm, if px1,then all m columns choosen to be same
%          Notes: nby{i,j} <= na{i}, strictly proper
%     solver_config: {solverName, lambdaList, timeSegments}
%          - solverName: string, "gl1", "gil1", "gsbl", "l1", "il1"; by default "gil1"
%          - lambdaList: double OR [double, double, double] <-> "grp_WE, grp_WS, E"
%          - timeSegments: int, OR, Kx2 vector:
%               int   : divide equally into K segements
%               vector: each row represents one segment;
%                       column 1 is the starting index,
%                       column 2 is the end index;
%
% Output:
%    narx_etime: the cost of cpu time in computation
%


% Globals for solver settings
global glb_optimizer
global MAX_ITER ABSTOL
config.name = glb_optimizer;
config.max_iter = MAX_ITER;
config.abstol = ABSTOL;


% Time Benchmark
arx_timer = tic;

% IO data dimensions
p = size(net.By, 1);
m = size(net.Bu, 2);

% Unpack datasets
if iscell(morders)
    na = morders{1};
    nby = morders{2};
    nbu = morders{3};
else  % 3x1 integer vector
    na = ones(p,1) * morders(1);
    nby = ones(p,p) * morders(2);
    nbu = ones(p,1) * morders(3);
end

% Standardize input dataset and settings
[na_sr, na_sc] = size(na);
[nby_sr, nby_sc] = size(nby);
[nbu_sr, nbu_sc] = size(nbu);
if na_sr ~= nby_sr || nby_sr ~= nby_sc || na_sc ~= 1
    error('The settings of morders is incorrect!')
end
if isempty(nbu) && (nbu_sr == 0 | nbu_sr == 1)
    nbu = zeros(p,0);
end
if na_sr == 1
    na = ones(p,1)*na;
    nby = ones(p,p)*nby;
    if ~isempty(nbu)
        nbu = ones(p,m)*nbu;
    end
end

if ~iscell(rawdata)
    % re-assemble the dataset according to "segments"
    timeseg = solver_config{3};
    if isscalar(timeseg)
        tEnd = size(rawdata.y, 1);
        timeseg = [linspace(1, round((tEnd-1)/timeseg)*timeseg+1, ...
                            timeseg+1), tEnd]';
        timeseg(end-1) = [];
        dataSetIndex = [timeseg(1:end-1), timeseg(2:end)];
    else
        dataSetIndex = timeseg;  % Kx2 vector, K segements
    end
    Cexp = size(dataSetIndex,1);
    data = cell(Cexp, 1);
    for i = 1:Cexp
        data{i} = rawdata(dataSetIndex(i,1):dataSetIndex(i,2));
    end
    clear('rawdata');
else
    error(['Boolean nonlinear dynamic network reconstruction has NOT ' ...
           'supported multiple experiments!'])
end
tM = size(data{1}.y, 1); % the length of time series



% Declare output variales for parfor
w_net = cell(p,1);
nBlockSizeList = cell(p,1);
pred_err = zeros(p,1);

parfor i = 1:p    % parfor
% -----------------------------------------------------
% Considering the i-th output y_i

% === Constructing the Regression form ===
if m   % having inputs
    orders = [na(i), nby(i,:), nbu(i,:)];
    nBlock = length(orders)-1;   % for N in 1:N
    nBlockSize = getBlockSizes(nBlock, [p,m,i], na(i), nby(i,:), nbu(i,:));
    blockIndexList = [0 cumsum(nBlockSize)];
    % the indices of each block is NblkIndexList(n)+1:NblkIndexList(n+1)
else   % no inputs
    orders = [na(i), nby(i,:)];
    nBlock = length(orders)-1;   % for N in 1:N
    nBlockSize = getBlockSizes(nBlock, [p,i], na(i), nby(i,:));
    blockIndexList = [0 cumsum(nBlockSize)];
end

% construct the essential Ac and phi matrices
cellA = {}; by = []; % by = Aw
for c = 1:Cexp
    % unpack data
    y = data{c}.y;
    u = data{c}.u;
    t1 = max(orders) + 1; % the index of t1
    % build A_i^c
    Ac = [];
    for t = t1:1:tM
        if m   % having inputs
            phi_t = buildPhiOI(t, [p,m,i], y,u, na(i),nby(i,:),nbu(i,:));
        else   % no inputs
            phi_t = buildPhiO(t, [p,i], y, na(i),nby(i,:));
        end
        Ac = [Ac; phi_t];
    end
    cellA = [cellA, Ac];
    by = [by; y(t1:tM,i)];
end

% construct the A matrix (by = A*w)
A = [];
for n = 1:nBlock
    eachBlock = blockIndexList(n)+1:blockIndexList(n+1);
    An = [];
    for c = 1:Cexp
        An = blkdiag(An, cellA{c}(:,eachBlock));
    end
    A = [A An];
end


% === Optimization ===
% wi = A\by;
% wiLasso = lasso(A, by);
% wi = wiLasso(:,85);

if length(solver_config{2}) == 1 || ...
   isscalar(solver_config{3}) && solver_config{3} == 1
%     warning('The lambda for small groups is set to ZERO!')
    lambda_WE = 0;
    lambda_WS = solver_config{2};
    lambda_E = Inf;
elseif length(solver_config{2}) == 2
    lambda_WE = solver_config{2}(1);
    lambda_WS = solver_config{2}(2);
    lambda_E = Inf;
else
    % Note: even if setting three diff lambdas, the specific sparsity
    % algorithm will only choose the one or two among them.
    lambda_WE = solver_config{2}(1); % small-group sparsity; used for comparing
    lambda_WS = solver_config{2}(2); % big-group sparsity; only use this one!
    lambda_E = solver_config{2}(3); % element sparsity; used for comparing
end

switch solver_config{1}
    case 'gl1'
        [w_net{i}, cod_k(i), rmse(i)] = GroupL1Solver_(by, A, ...
                                  {[nBlock, Cexp], nBlockSize},...
                                  [lambda_WE, lambda_WS, lambda_E],...
                                  config);
        % solver_name = 'Group LASSO';
    case 'gil1'  % support
        [w_net{i}, pred_err(i)] = NGroupIL1Solver_(by, A, ...
                                   {[nBlock, Cexp], nBlockSize},...
                                   [lambda_WE, lambda_WS, lambda_E],...
                                   config);
        % solver_name = 'Group Iterative Reweighted l1 Method';
    case 'gil2'
        [w_net{i}, cod_k(i), rmse(i)] = GroupIL2Solver_(by, A, ...
                                 {[nBlock, Cexp], nBlockSize},...
                                 [lambda_WE, lambda_WS, lambda_E]);
        % solver_name = 'Group Iterative Reweighted l2 Method';
    case 'gsbl'
        [w_net{i}, cod_k(i), rmse(i)] = GroupSBLSolver_(by, A, ...
                                 {[nBlock, Cexp], nBlockSize},...
                                 [lambda_WE, lambda_WS, lambda_E]);
        % solver_name = 'Sparse Bayesian Learning (group sparsity)';
    case 'l1'
        [w_net{i}, cod_k(i), rmse(i)] = L1Solver_(by, A, ...
                                 {[nBlock, Cexp], nBlockSize},...
                                 [lambda_WE, lambda_WS, lambda_E]);
        % solver_name = 'LASSO';
    case 'il1'
        [w_net{i}, cod_k(i), rmse(i)] = IL1Solver_(by, A, ...
                                 {[nBlock, Cexp], nBlockSize},...
                                 [lambda_WE, lambda_WS, lambda_E]);
        % solver_name = 'Iterative Reweighted l1 Method';
    case 'il2'
        [w_net{i}, cod_k(i), rmse(i)] = IL2Solver_(by, A, ...
                                 {[nBlock, Cexp], nBlockSize},...
                                 [lambda_WE, lambda_WS, lambda_E]);
        % solver_name = 'Iterative Reweighted l2 Method';
end

% === Save solutions over parfor loops ===
% w_net{i}, pred_err(i)
nBlockSizeList{i} = nBlockSize;

% -----------------------------------------------------
end  % END of i-loop

% === reconstruct parameters from w and save in idnet obj
arxA = cell(1,p); arxBy = cell(p,p); arxBu = cell(p,m);
for i = 1:p
    nBlock = size(nBlockSizeList{i},2);
    [arxA(i), arxBy(i,:), arxBu(i,:)] = w2AByBu(w_net{i}, ...
        [p,m,i], [Cexp, nBlock], nBlockSizeList{i});
end
net.A = arxA; net.By = arxBy; net.Bu = arxBu;

% === writing other class members
net.Ts = data{1}.Ts;
net.BoolNet = net.bool('NL');
net.PredError = pred_err;
net.Report.Solver = SolverNameMap(solver_config{1});
net.Report.Optimizer = glb_optimizer;


% === print out the CPU time used
narx_etime = toc(arx_timer);
% fprintf('The elapsed time for arx.m: %f seconds \n', narx_etime);

end % END of Function







% ====================================
% Local functions
% ====================================
function phiT = buildPhiOI(t, pmi, y, u, nai, nbyi, nbui)
% if nby(i,j) == const, Forall j; so nbu(i,j)
% there is a faster matrix manipulation.
    p = pmi(1); m = pmi(2); i = pmi(3);
    phiT = [];
    for j = 1:p
        if j == i
            phiT = [phiT -flipud(y(t-nai:t-1, j)')];
        else
            phiT = [phiT flipud(y(t-nbyi(j):t-1, j)')];
        end
    end
    for j = 1:m
        phiT = [phiT flipud(u(t-nbui(j):t-1, j)')];
    end
end

function phiT = buildPhiO(t, pmi, y, nai, nbyi)
    p = pmi(1); i = pmi(2);
    phiT = [];
    for j = 1:p
        if j == i
            phiT = [phiT -flipud(y(t-nai:t-1, j)')];
        else
            phiT = [phiT flipud(y(t-nbyi(j):t-1, j)')];
        end
    end
end

function [rowA, rowBy, rowBu] = w2AByBu(w, pmi, CN, nBlockSize)
    p = pmi(1); m = pmi(2); i = pmi(3);
    C = CN(1); N = CN(2);
    IndexMapWS2 = indexMap([C,N], 'ws2', nBlockSize);

    % declare variables to save data, then save to object net
    rowA = cell(1,1); rowBy = cell(1,p); rowBu = cell(1,m);

    % re-assemble: w -> A, By, Bu
    for k = 1:N
        index_ws2 = k;
        n_mgrp = IndexMapWS2{index_ws2};  % indices of multiple experiments
        n_mgrp = reshape(n_mgrp, [], C);

        wList = {}; zeroCell = {};
        for c = 1:C
            n = n_mgrp(:,c);
            if k == i
                wList = [wList [1 w(n)']];
                zeroCell = [zeroCell 0];
            else
                wList = [wList w(n)'];
            end
        end
        if 1 <= k && k <= p
            rowBy{k} = wList;
            if k == i
                rowBy{k} = zeroCell;
                rowA{1} = wList;
            end
        else  % if having inputs
            rowBu{k-p} = wList;
        end
    end
end

function indexList = getBlockSizes(N, pmi, na, nby, nbu)
    if length(pmi) == 3
        p = pmi(1); m = pmi(2); i = pmi(3);
        indexList = [];
        for k = 1:N
            if k == i
                indexList = [indexList na];
            elseif k <= p
                indexList = [indexList nby(k)];
            elseif k > p
                indexList = [indexList nbu(k-p)];
            end
        end
    elseif length(pmi) == 2
        p = pmi(1); i = pmi(2);
        indexList = [];
        for k = 1:N
            if k == i
                indexList = [indexList na];
            else
                indexList = [indexList nby(k)];
            end
        end
    end
end


% index maps from w^S1, w^S2 to w
function wIndex = indexMap(CN, varname, blkSize)
    C = CN(1); N = CN(2);
    NS1 = C*N; NS2 = N;
    if strcmp(varname, 'ws1')  % small groups; no longer useful!
        for k = 1:NS1
            wIndex{k} = C * sum(blkSize(1:ceil(k/C)-1)) * ones(1,blkSize(ceil(k/C))) ...
                + (rem(k-1,C) * blkSize(ceil(k/C))+1 : (rem(k-1,C)+1)*blkSize(ceil(k/C)));
        end
    elseif strcmp(varname, 'ws2')
        for k = 1:NS2
            wIndex{k} = C * sum(blkSize(1:k-1)) * ones(1, C*blkSize(k)) ...
                + (1:C*blkSize(k));
        end
    end
end