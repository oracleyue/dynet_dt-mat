function [w_est, sigma2_est, gamma_est, solver_status] = ...
    groupSBL(Phi, y, mGroupTable, lambda, optimizerOpt)
% GROUPSBL uses Sparse Bayesian Learning for sparsity solutions to y = Phi*w.
%
% INPUT:
%   Phi                 : double, Ny x Nw vector
%   y                   : double, Ny x 1 vector
%   mGroupindex         : integer, N x 2 vector (N      : number of groups)
%   lambda              : double
%   optimizerOpt        : struct, to control optimization process
%       - name (not used)
%       - MAXITER: maximal iterations (default: 1e4)
%       - ABSTOL: absolute tolerance (default: 1e-8)
%       - RELTOL (not used)
%       - PruneGamma: threshold to prune gamma (default: 1e-4)
%
% OUTPUT:
%   w_est               : parameter estimation
%   sigma2_est          : covariance matrix of noise: sigma^2 I
%   gamma_est           : hyperparameters, the covariance of Gaussian prior
%
%
% REFERENCES:
% [1] Wipf, D. P., & Rao, B. D. (2004). Sparse Bayesian learning for basis selection. Signal Processing, IEEE Transactions on, 52(8), 2153-2164.
% [2] Yue, Z., Pan, W., Thunberg, J., Ljung, L., & Goncalves, J. (2017). Linear Dynamic Network Reconstruction from Heterogeneous Datasets. In Preprints of the 20th World Congress, IFAC (pp. 11075–11080). Toulouse, France.
%
% ACKNOWLEDGEMENT:
%   The code "MSBL.m" written by Zhilin Zhang.

% Copyright (c) 2015-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 16 Aug 2017



%
% Dimension of the Problem
%
[Ny, Nw]   = size(Phi);          % Ny,Nw: const
[Ngrp, ~]  = size(mGroupTable);  % Ngrp: const; mVar: matrix
vGroupSize = mGroupTable(:,2) - mGroupTable(:,1) + 1;    % vVar: vector
Mw         = sum(vGroupSize);    % we update Mw when removing zeros from "w"
% check dimension compatibility
if Ny ~= length(y);
    error('The dimensions of y and Phi fail to match!');
elseif Nw ~= Mw
    error(['The dimensions of Phi and mGroupTable (group index table)' ...
           'fail to match!']);
end


%
% Default Control Parameters
%
if nargin == 4
    PRUNE_GAMMA = 1e-4;    % threshold for prunning small hyperparameters gamma_i
    EPSILON     = 1e-8;    % threshold for stopping iteration.
    MAX_ITER    = 2000;    % maximum iterations
    LearnLambda = 1;       % control update sigma^2 or not
elseif nargin == 5
    PRUNE_GAMMA = optimizerOpt.PruneGamma;
    EPSILON     = optimizerOpt.ABSTOL;
    MAX_ITER    = optimizerOpt.MAXITER;
    LearnLambda = optimizerOpt.LearnLambda;
else
    error('Too many or not enough arguments!')
end
PRINT       = 0;       % for debugging

if (PRINT) fprintf('\nRunning SBL ...\n'); end


%
% Initializations
%
Mgrp           = Ngrp;            % number of groups after pruning
sigma2         = lambda;          % sigma^2
gamma          = ones(Mgrp,1);    % hyperparameter
gammaIndexList = (1:Mgrp)';       % indexex of non-zero hyperparameter
mu             = zeros(Mw,1);     % ML mean of w
count          = 0;               % iteration count


%
% Algorithm Implementations (EM)
%
while 1

    %
    % Prune weights as their hyperparameters go to zero
    %
    if min(gamma) < PRUNE_GAMMA
        index = find(gamma > PRUNE_GAMMA);
        Mgrp = length(index);
        gammaIndexList = gammaIndexList(index);
        vIndexList = ...
            groupTable2indexList(updateGroupTable(mGroupTable), index);
        gamma = gamma(index);  % use elements larger than MIN_GAMMA
        Phi = Phi(:,vIndexList);  % corresponding columns in Phi
        mGroupTable = mGroupTable(index, :);
        vGroupSize = vGroupSize(index);
        Mw = sum(vGroupSize);
    end

    mu_old =mu;
    % Gamma = blkdiag(gamma(1)I_1, ..., gamma(m)I_m)
    GammaDiag = expandVariable(gamma, vGroupSize);
    Gamma = diag(GammaDiag);
    G = diag(sqrt(GammaDiag));

    %
    % Update the mean $\mu$ and covariance $\sigma^2$
    %
    % refer to Eq.(18),(19) in [1] but modified for the grouped version
    [U,S,V] = svd(Phi*G,'econ');
    diag_S = diag(S);
    Xi = G * V * diag((diag_S./(diag_S.^2 + sigma2 + 1e-16))) * U';
    mu = Xi * y;

    %
    % Update hyperparameters
    %
    % refer to Eq.(M-step in [2]); mu2 := mu.^2, SigmaWDiag := (\Sigma_w)_ii
    mu2 = mu.^2;
    % SigmaWDiag = real(GammaDiag - (sum(Xi'.*(Phi*Gamma)))');   % MSB.m version
    % SigmaWDiag = GammaDiag - sum(Xi*Phi*Gamma,2);   % Sigma_ii <- mean(Sigma(i,:))
    SigmaWDiag = GammaDiag - diag(Xi*Phi*Gamma);      % Sigma_ii

    mu2_group        = zeros(Mgrp,1);    % average of mu2 over each group
    SigmaWDiag_group = zeros(Mgrp,1);    % average of SigmaWDiag over each group
    mGroupTableLocal = updateGroupTable(mGroupTable);
    % taking average of mu2 and diag(sigma_w) over each group
    for indexGroup = 1:Mgrp
        vIndex = mGroupTableLocal(indexGroup,1):1:mGroupTableLocal(indexGroup,2);
        groupSize = vGroupSize(indexGroup);
        mu2_group(indexGroup) = sum(mu2(vIndex))/groupSize;
        SigmaWDiag_group(indexGroup) = sum(SigmaWDiag(vIndex))/groupSize;
    end
    gamma_old = gamma;
    gamma = mu2_group + SigmaWDiag_group;

    %
    % Update sigma^2, i.e. the learning rule for lambda
    %
    % Note: You can use it to estimate the lambda when SNR >= 20 dB. But when SNR <
    % 20 dB, you'd better use other methods to estimate the lambda, since it is not
    % robust in strongly noisy cases (in simulations, you can feed it with the true
    % noise variance, which can lead to a near-optimal performance)
    if LearnLambda == 1
        sigma2_old = sigma2;
        % Rule 1: in general
        % sigma2 = norm(y-Phi*mu, 2)^2 / Ny + ...
        %          sigma2_old * (Mw - sum(SigmaWDiag_group./gamma_old ...
        %                           .* vGroupSize)) / Ny;

        % Rule 2: force to zero
        sigma2 = norm(y-Phi*mu,2)^2 / (Ny - Mw + ...
                 sum(SigmaWDiag_group./gamma_old .* vGroupSize));
    end

    %
    % Stopping conditions
    %
    count = count + 1;
    if (PRINT) disp(['iters: ',num2str(count),'   num coeffs: ',num2str(Mw), ...
                     '   gamma change: ',num2str(max(abs(gamma - gamma_old))), ...
                     '   dim(mu) : ',num2str(length(mu))]); end

    if count >= MAX_ITER
        solver_status = 'stop on MAX_ITER'; break
    end
    if size(mu) == size(mu_old)
        diff_mu = max(max(abs(mu_old - mu)));
        if diff_mu < EPSILON
            solver_status = 'stop on \\mu'; break
        end
        % diff_gamma = max(max(abs(gamma_old - gamma)));
        % if diff_gamma < EPSILON
        %     solver_status = 'stop on \\gamma'; break
        % end
    end

end


%
% Expanding solutions
%

% noise variance
sigma2_est = sigma2;

% hyperparameters "gamma"
gamma_est = zeros(Ngrp,1);
gamma_est(gammaIndexList) = gamma;

% parameters "w"
vElementIndex = ...       % indexex of non-zero parameter
    groupTable2indexList(mGroupTable, 1:size(mGroupTable,1));
w_est = zeros(Nw,1);
w_est(vElementIndex,:) = mu;

if (PRINT) fprintf('solver status: %s\n', solver_status); end
if (PRINT) fprintf('\nFinish running.\n'); end

end % END:main


% ================================================================
% Local Functions
% ================================================================

% Retrieve elementary indexes given the group table.
% Info: a group index table is (Mgrp x 2), in each row the first element
% is the starting index of that group, and the second the ending index.
% e.g.
%     mGroupTable  [3   5;   ==>  vIndexList   [3 4 5 10 11 12 13]
%                   10 13]
function vIndexList = groupTable2indexList(mGroupTable, groupIndexSelected)
    vIndexList = [];
    for k = 1:length(groupIndexSelected)
        index = groupIndexSelected(k);
        vIndexList = [vIndexList ...
                      mGroupTable(index,1):1:mGroupTable(index,2)];
    end
end

% Update the original group indexes of w to the group indexes of parameters
% after removing zeros, e.g.
%                   [4   6;                                 [1  3;
%     mGroupTable    8  10;      ==>    mGroupTableLocal     4  6;
%                   12  15]                                  7  10]
%
function mGroupTableLocal = updateGroupTable(mGroupTable)
    mGroupTableLocal = zeros(size(mGroupTable,1),2);
    initIndex = 1;
    for index = 1:size(mGroupTable,1)
        mGroupTableLocal(index,1) = initIndex;
        mGroupTableLocal(index,2) = initIndex + ...
            mGroupTable(index,2) - mGroupTable(index,1);
        initIndex = mGroupTableLocal(index,2) + 1;
    end
end

% Expand the grouped version variable to the elementary verion,
% e.g. expand gamma (group version) to gamma (element version), i.e.
%   Gamma = blkdiag(gamma(1)*I_1, ..., gamma(Mgrp)*I_{Mgrp})
%   gamma (element version):= diag(Gamma)
function vVariableElement = expandVariable(vVariableGroup, vGroupSize)
    vVariableElement = [];
    for indexGroup = 1:length(vGroupSize)
        varBlock = vVariableGroup(indexGroup)*ones(vGroupSize(indexGroup),1);
        vVariableElement = [vVariableElement; varBlock];
    end
end
