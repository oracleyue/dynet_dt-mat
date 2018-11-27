addpath('../solvers/')

renew_data = 0;
if renew_data seed = rng; else rng(seed); end

N = 100;   % num of observations
dimG = [2 4 3 2 3 3 2 3 2];
numZeroGroup = 4;
M = sum(dimG);   % dim of parameters
snr = 40; % unit: dB   10log(varSignal/varNoise)

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

algorithm = 'gspmc';
mcpara.iter_max = 1e4;
switch algorithm
  case 'gspmc'
    indexLists{1}(1) = length(dimG);
    indexLists{1}(2) = 1;
    indexLists{2} = dimG;
    solverOptions = solverInit('gspmc', .1);
    [w_est, gtopol, gprob, nvar] = ...
        GroupSpMCSolver_(t, Phi, indexLists, solverOptions.optimizerOpt);
    disp([w w_est])
    disp([spGroup gtopol gprob])

  case 'gsbl'
    mGroupTable = [cumsum([1 dimG(1:end-1)]); cumsum(dimG)];
    [w_est, sigma2_est, gamma_est, SBL_SolverStatus] = ...
        groupSBL(Phi, t, mGroupTable', lambda);
    fprintf('Ground truth w and w_est:\n');
    disp([w w_est])



    disp([w w_est])

end