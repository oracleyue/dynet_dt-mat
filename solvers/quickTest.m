renew_data = 1;
if renew_data seed = rng; else rng(seed); end


N = 1000;   % num of observations
dimG = [2 4 3 2 3 3 2 3 2];
numZeroGroup = 4;
M = sum(dimG);   % dim of parameters
snr = 10;   % unit: dB   10log(varSignal/varNoise)

Phi = randn(N, M);
w = randn(M, 1);
indices = randi(length(dimG), numZeroGroup, 1);
sgt = ones(length(dimG),1);
for k = 1:numZeroGroup
    index = indices(k);
    if index == 1
        lIndex = 1;
    else
        lIndex = sum(dimG(1:index-1))+1;
    end
    rIndex = sum(dimG(1:index));
    w(lIndex:rIndex) = 0;
    sgt(index) = 0;
end
y = Phi * w;
varY = std(y)^2;
varNoise = varY / (10^(snr/10));
noise = randn(N,1) * sqrt(varNoise);
t = y+noise;

% select time interval
ts_range = 1:100;
t = t(ts_range);
Phi = Phi(ts_range,:);

sTimer = tic;
algorithm = 'g';
hypara.R0 = 1e-1;
mcpara.iter_max = 4e3;
mcpara.ENABLE_DEBUG = 1;
switch algorithm
  case {'EGHB', 'eghb'}
    [s, sG, R, gamma] = fEGSparseHBMC(t, Phi, dimG, 'mcpara', mcpara);
    fprintf('Ground truth w and prob of element & group sparsity s:\n');
    disp([w s [sG; NaN*ones(length(s)-length(sG),1)]])

  case {'EGUB', 'egub'}
    [s, sG, R, gamma] = fEGSparseUBMC(t, Phi, dimG, 'mcpara', mcpara);
    fprintf('Ground truth w and prob of element & group sparsity s:\n');
    disp([w s [sG; NaN*ones(length(s)-length(sG),1)]])

  case {'EGC', 'egc'}
    [s, sG, R, gamma] = fEGCSparseMC(t, Phi, dimG, 'mcpara', mcpara);
    fprintf('Ground truth w and prob of element & group sparsity s:\n');
    disp([w s [sG; NaN*ones(length(s)-length(sG),1)]])

  case {'EG', 'eg'}
    [s, sG, R, gamma] = fEGSparseMC(t, Phi, dimG, 'mcpara', mcpara);
    fprintf('Ground truth w and prob of element & group sparsity s:\n');
    disp([w s [sG; NaN*ones(length(s)-length(sG),1)]])

  case {'G', 'g'}
    [s, R, gamma] = groupSpMC(t, Phi, dimG, ...
                              'mcpara', mcpara, 'hypara', hypara);
    fprintf('Ground truth w and prob of group sparsity s:\n');
    disp([sgt s])

  case {'E', 'e'}
    [sG, R, gamma] = fESparseMC(t, Phi, 'mcpara', mcpara);
    fprintf('Ground truth w and prob of element sparsity s:\n');
    disp([w s])

  case 'gsbl'
    mGroupTable = [cumsum([1 dimG(1:end-1)]); cumsum(dimG)];
    [w_est, sigma2_est, gamma_est, SBL_SolverStatus] = ...
        groupSBL(Phi, t, mGroupTable', lambda);
    fprintf('Ground truth w and w_est:\n');
    disp([w w_est])

  case 'sbl'
    [w_est, sigma2_est, gamma_est, SBL_SolverStatus] = ...
        MSBL(Phi, t, lambda, 1);
    fprintf('Ground truth w and w_est:\n');
    disp([w w_est])

end
toc(sTimer)