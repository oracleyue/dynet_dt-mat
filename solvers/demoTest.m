renew_data = 1;
if renew_data seed = rng; else rng(seed); end

% sim data
N = 1000;   % num of observations
dimG = [2 4 3 2 3 3 2 3 2];
numZeroGroup = 4;
M = sum(dimG);   % dim of parameters
Mg = length(dimG);
snr = 10; % unit: dB   10log(varSignal/varNoise)

Phi = randn(N, M);
w = randn(M, 1);
indices = randi(length(dimG), numZeroGroup, 1);
spGroup = ones(length(dimG),1);
for k = 1:numZeroGroup
    index = indices(k);
    if index == 1
        lIndex = 1;
    else
        lIndex = sum(dimG(1:index-1))+1;
    end
    rIndex = sum(dimG(1:index));
    w(lIndex:rIndex) = 0;
    spGroup(index) = 0;
end
y = Phi * w;
varY = std(y)^2;
varNoise = varY / (10^(snr/10));
noise = randn(N,1) * sqrt(varNoise);
t = y+noise;

% estimation
sTimer = tic;
algorithm = 'gspmc';
switch algorithm
  case 'gspmc'
    hypara.a = 1e-4;
    hypara.b = 1e-4;
    hypara.c = 1e-4;
    hypara.d = 1e-4;
    hypara.pB = .5;
    hypara.R0 = 1.6;
    hypara.gamma0 = 1; % ones(Mg,1);
    mcpara.iter_max = 4e3;  % max number of Gibbs iterations
    mcpara.var_gamma = .01;
    mcpara.var_R = 1e-2;
    mcpara.prunePhi = 0;
    mcpara.MH_scheme = 1;  % 1 or 2; schemes of proposal dist of s
    mcpara.Pred_scheme = 'mean';
    mcpara.ENABLE_DEBUG = 1;

    [sG, R, gamma] = groupSpMC(t, Phi, dimG, ...
                              'hypara', hypara, 'mcpara', mcpara);
    disp([spGroup sG])
    disp([varNoise, R])

  case 'gsbl'
    mGroupTable = [cumsum([1 dimG(1:end-1)]); cumsum(dimG)];
    [w_est, sigma2_est, gamma_est, SBL_SolverStatus] = ...
        groupSBL(Phi, t, mGroupTable', lambda);
    fprintf('Ground truth w and w_est:\n');
    disp([w w_est])

end
toc(sTimer)