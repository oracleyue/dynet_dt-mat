function [w_est, lambda_est, gamma_est, solver_status] = ...
    classicSBL(Phi, y, lambda, optimizerOpt)
% ELEMENTSBL uses Sparse Bayesian Learning for sparsity solutions to y = Phi*w.
%
% INPUT:
%   Phi                 : double, Ny x Nw vector
%   y                   : double, Ny x 1 vector
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
%   lambda_est          : covariance matrix of noise: sigma^2 I
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
% Last update on 23 Jul 2018


% Dimension of the Problem
[Ny Nw] = size(Phi);
if Ny ~= length(y);
    error('The dimensions of y and Phi fail to match!');
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
sigma2         = lambda;          % sigma^2
gamma          = ones(Nw,1);      % hyperparameter
keep_list      = (1:Nw)';         % indexex of non-zero hyperparameter
mu             = zeros(Nw,1);     % ML mean of w
count          = 0;               % iteration count
m = length(keep_list);


%
% EM Algorithm
%
while 1

    % *** Prune weights as their hyperparameters go to zero ***
    if min(gamma) < PRUNE_GAMMA
        index = find(gamma > PRUNE_GAMMA);
        gamma = gamma(index);  % use all the elements larger than MIN_GAMMA to form new 'gamma'
        Phi = Phi(:,index);    % corresponding columns in Phi
        keep_list = keep_list(index);
        m = length(gamma);
    end

    mu_old =mu;
    Gamma = diag(gamma);
    G = diag(sqrt(gamma));

    % ****** estimate the solution matrix *****
    [U,S,V] = svd(Phi*G,'econ');

    [d1,d2] = size(S);
    if (d1 > 1)     diag_S = diag(S);
    else            diag_S = S(1);      end

    Xi = G * V * diag((diag_S./(diag_S.^2 + lambda + 1e-16))) * U';
    mu = Xi * y;

    % *** Update hyperparameters, i.e. Eq(18) in the reference ***
    gamma_old = gamma;
    mu2_bar = sum(abs(mu).^2,2)/L;

    Sigma_w_diag = real( gamma - (sum(Xi'.*(Phi*Gamma)))');
    gamma = mu2_bar + Sigma_w_diag;

    % ***** the lambda learning rule *****
    % You can use it to estimate the lambda when SNR >= 20 dB. But when SNR < 20 dB,
    % you'd better use other methods to estimate the lambda, since it is not robust
    % in strongly noisy cases (in simulations, you can feed it with the
    % true noise variance, which can lead to a near-optimal performance)
    if LearnLambda == 1
        lambda = (norm(y - Phi * mu,'fro')^2/L)/(Ny-m + sum(Sigma_w_diag./gamma_old));
    end

    % *** Check stopping conditions, etc. ***
    count = count + 1;
    if (PRINT) disp(['iters: ',num2str(count),'   num coeffs: ',num2str(m), ...
            '   gamma change: ',num2str(max(abs(gamma - gamma_old)))]); end
    if (count >= MAX_ITERS)
        solver_status = 'stop on MAX_ITER'; break
    end

    if (size(mu) == size(mu_old))
        dmu = max(max(abs(mu_old - mu)));
        if dmu < EPSILON
            solver_status = 'stop on \\mu'; break
        end
        % diff_gamma = max(max(abs(gamma_old - gamma)));
        % if diff_gamma < EPSILON
        %     solver_status = 'stop on \\gamma'; break
        % end
    end

end


%
% Expand solutions
%

% hyperparameters
gamma_ind = sort(keep_list);
gamma_est = zeros(Nw,1);
gamma_est(keep_list,1) = gamma;

% noise variance
lambda_est = lambda;

% expand the final solution
w_est = zeros(Nw,1);
w_est(keep_list,:) = mu;

if (PRINT)
    fprintf('solver status: %s\n', solver_status);
    fprintf('\nFinish running ...\n');
end

return % END:main
