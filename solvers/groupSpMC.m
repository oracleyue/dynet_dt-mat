function [sG, R, gamma] = groupSpMC(t, Phi, dimG, varargin)
% GROUPSPMC to infer group-wise sparse structure of w
% from t = Phi*w + N(0,R*I).
% INFO: This is a modified version (name, interface) of "fGSparseMC.,".
%
% This function considers the simple case (single variate regression), using
% scalar R and gamma, and hence scalar hyperparameters (a,b,c,d). Each
% element of t is an independent observation t_i = Phi(i,:)*w + N(0,R).
%
% INPUT:
%   t       : Nx1 vector; dependent variables
%   Phi     : NxM matrix; regressors
%   dimG    : Ngx1 vector: group sizes (sum(dimG) == M)
%   hypara  : (optional) struct; hyper parameters
%     - a,b   : R ~ invGamma(a,b); suggested: 1e-4 (both)
%     - c,d   : gamma ~ invGamma(c,d); suggested: 1e-4 (both)
%     - pB : Bernoulli prior for s_i ~ B(1,pB); default: .5
%     (if given, use them to initialize R and gamma)
%     - R     : scalar; noise variacne
%     - gamma : scalar (default) OR Mx1 vector; variance of prior of w
%   mcpara  : (optional) struct; tuning MCMC performance
%     - iter_max : max number of MC iterations; default: 1e4
%     - var_gamma   : random walk for gamma; default: .01
%     - var_R   : random walk for R; default: .01
%     - prunePhi : threshold to prune Phi based on s/sG when nonzero
%     - MH_scheme : {1, 2}; specify the proposal distribution for s
%     - Pred_scheme : {'mean', 'prob'}; choose s by using the mean of
%     samples or choose by the appearance frequencies among samples.
%     - ENABLE_DEBUG : 1 to print  debug information
%
% OUTPUT:
%   s      :  grouply sparse structure of the parameter vector
%   R      :  noise variance
%   gamma  :  hyperparameter; w ~ N(0,gamma*I) or N(0,diag(gamma))
%
% Examples:
%
%   Demo 1: the simplest usage
%   [s, R, gamma] = fSpMCLite(t, Phi)
%
%   Demo 2: adjust hyperprior of R and gamma
%   hypara.a = .01; hypara.b = .01; hypara.c = .01; hypara.d = .01;
%   [s, R, gamma] = fSpMCLite(t, Phi, 'hypara', hypara);
%
%   Demo 3: prefer sparser structure (set pB smaller)
%   hypara.pB = .02;
%   [s, R, gamma] = fSpMCLite(t, Phi, 'hypara', hypara);
%
%   Demo 4: use diagonal Gamma instead of scalar
%   hypara.gamma = ones(size(Phi,2),1);
%   [s, R, gamma] = fSpMCLite(t, Phi, 'hypara', hypara);
%
%   Demo 5: tune for better MCMC performance
%   mcpara.iter_max = 1e5;
%   mcpara.Pred_scheme = 'prob';
%   mcpara.var_gamma = .1;
%   [s, R, gamma] = fSpMCLite(t, Phi, 'mcpara', mcpara);
%   OR, simply
%   [s, R, gamma] = fSpMCLite(t, Phi, mcpara);

% Copyright (c) 2018, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 15 Jun 2018


% Parsing function arguments
argParser = inputParser;
argParser.FunctionName = 'fSpMCLite';

hyparaDefault.a = 1e-4;
hyparaDefault.b = 1e-4;
hyparaDefault.c = 1e-4;
hyparaDefault.d = 1e-4;
hyparaDefault.pB = .5;
hyparaDefault.R0 = 1;
hyparaDefault.gamma0 = 1;
mcparaDefault.iter_max = 1e4;  % max number of Gibbs iterations
mcparaDefault.var_gamma = .01;
mcparaDefault.var_R = .01;
mcparaDefault.prunePhi = 0;   % no prune action
mcparaDefault.MH_scheme = 1;  % 1 or 2; schemes of proposal dist of s
mcparaDefault.Pred_scheme = 'mean';
mcparaDefault.ENABLE_DEBUG = 0;

addRequired(argParser, 't');
addRequired(argParser, 'Phi');
addRequired(argParser, 'dimG');
addParameter(argParser, 'mcpara', mcparaDefault);
addParameter(argParser, 'hypara', hyparaDefault);
parse(argParser, t, Phi, dimG, varargin{:});
t = argParser.Results.t;
Phi = argParser.Results.Phi;
dimG = argParser.Results.dimG;
hypara = argParser.Results.hypara;
mcpara = argParser.Results.mcpara;

% tuning for better MCMC performance
if isfield(mcpara, 'iter_max')
    iter_max = mcpara.iter_max;  % max number of gibbs iterations
else
    iter_max = mcparaDefault.iter_max;
end
if isfield(mcpara, 'var_gamma')
    var_gamma = mcpara.var_gamma;
else
    var_gamma = mcparaDefault.var_gamma;
end
if isfield(mcpara, 'var_R')
    var_R = mcpara.var_R;
else
    var_R = mcparaDefault.var_R;
end
if isfield(mcpara, 'prunePhi')
    prunePhi = mcpara.prunePhi;
else
    prunePhi = mcparaDefault.prunePhi;
end
if isfield(mcpara, 'MH_scheme')
    MH_scheme = mcpara.MH_scheme;
else
    MH_scheme = mcparaDefault.MH_scheme;
end
if isfield(mcpara, 'Pred_scheme')
    Pred_scheme = mcpara.Pred_scheme;
else
    Pred_scheme = mcparaDefault.Pred_scheme;
end
if isfield(mcpara, 'ENABLE_DEBUG')
    ENABLE_DEBUG = mcpara.ENABLE_DEBUG;
else
    ENABLE_DEBUG = mcparaDefault.ENABLE_DEBUG;
end

% hyperparameters, initialization of random variables, etc.
if isfield(hypara, 'a') && isfield(hypara, 'b') && ...
        isfield(hypara, 'c') && isfield(hypara, 'd')
    a = hypara.a; b = hypara.b; c = hypara.c; d = hypara.d;
else
    a = hyparaDefault.a; b = hyparaDefault.b;
    c = hyparaDefault.c; d = hyparaDefault.d;
end
if isfield(hypara, 'pB')
    pB = hypara.pB;
else
    pB = hyparaDefault.pB;
end
if isfield(hypara, 'R0')
    R0 = hypara.R0;
else
    R0 = hyparaDefault.R0;
end
if isfield(hypara, 'gamma0')
    gamma0 = hypara.gamma0;
else
    gamma0 = hyparaDefault.gamma0;
end

% check dimension compatiability
N = length(t);
[Nphi, M] = size(Phi);
if N ~= Nphi
    error('The dimensions of t and Phi fail to match.')
end

if isfield(hypara, 'R')
    if size(R,1) ~= 1
        error('R needs to be a positive real scalar.')
    end
end

if isfield(hypara, 'gamma')
    if size(gamma,1) ~= 1 && size(gamma,1) ~= M
        error(['gamma can either be a positive real scalar or a real ' ...
               'positive vector with the dimension of unknown parameters.'])
    end
end


% Gibbs sampler + Metropolis-Hastings

% debugging start
if ENABLE_DEBUG
    fprintf('MCMC sampling in progress:\n');
end


% initialization
Ng = length(dimG);
sG_k = zeros(Ng,1);
gamma_k = gamma0;
R_k = R0;

% cumulate samples in MCMC processes
sampl_sG = sG_k;
sampl_gamma = gamma_k;
sampl_R = R_k;

if ENABLE_DEBUG
    MHAcceptRatio = ones(3,1);
    MHProbVal = zeros(3,1);
end

% Gibbs loop
for iter = 1:iter_max

    % draw s
    switch MH_scheme
      case 1
        % proposal distrition
        i_rand = randi(Ng);
        sG_hat = sG_k;  % \hat{s}
        sG_hat(i_rand) = 1 - sG_hat(i_rand);
        % s_k = groupLift(sG_k, dimG);
        % s_hat = groupLift(sG_hat, dimG);

        % M-H algorithm
        U = rand(1);
        pval_sk = evalProb(sG_k, gamma_k, R_k);
        pval_shat = evalProb(sG_hat, gamma_k, R_k);
        acceptRatio_s = min(1,  exp(pval_shat-pval_sk));
        if U <= acceptRatio_s
            sG_k = sG_hat;
        end
        sampl_sG = [sampl_sG sG_k];

      case 2
        % to be done
    end

    % draw gamma
    if length(gamma0) == 1  % use one sacalar gamma for every element of w
        gamma_hat = gamma_k + randn(1,1)*sqrt(var_gamma);
    else  % use diagonal Gamma for prior of w
        gamma_hat = gamma_k + randn(Ng,1)*sqrt(var_gamma);
    end
    gamma_hat = abs(gamma_hat);  % mirror negative values
    U = rand(1);
    pval_ghat = evalProb(sG_k, gamma_hat, R_k);
    pval_gk = evalProb(sG_k, gamma_k, R_k);
    acceptRatio_g = min(1, exp(pval_ghat-pval_gk));
    if U <= acceptRatio_g
        gamma_k = gamma_hat;
    end
    sampl_gamma = [sampl_gamma gamma_k];

    % draw R (scalar R)
    R_hat = R_k + randn(1)*sqrt(var_R);
    R_hat = abs(R_hat);   % mirror negative values
    U = rand(1);
    pval_Rhat = evalProb(sG_k, gamma_k, R_hat);
    pval_Rk = evalProb(sG_k, gamma_k, R_k);
    acceptRatio_r = min(1, exp(pval_Rhat-pval_Rk));
    if U <= acceptRatio_r
        R_k = R_hat;
    end
    sampl_R = [sampl_R R_k];

    if ENABLE_DEBUG
        MHAcceptRatio = [MHAcceptRatio ...
                         [acceptRatio_s; acceptRatio_g; acceptRatio_r]];
        MHProbVal = [MHProbVal ...
                     [pval_sk; pval_gk; pval_Rk]];
    end

   % debugging output
   if ENABLE_DEBUG && ~rem(iter, round(iter_max/50))
       fprintf('... iter %s:\n', num2str(iter));
       fprintf('    accept ratios: %d, %d, %d\n', ...
               acceptRatio_s, acceptRatio_g, acceptRatio_r);
       fprintf('    zero norm of sG: %d\n', sum(sG_k));
   end
end

range = round(iter_max/20):size(sampl_sG,2);  % discard first 20% samples
if strcmp(Pred_scheme, 'mean')
    sG = mean(sampl_sG(:,range), 2);
elseif strcmp(Pred_scheme, 'prob')
    sG = findMaxFreq(sampl_sG(:,range));
end
% s = groupLift(sG, dimG);
% return sG instead of s
if isscalar(gamma0)
    gamma = mean(sampl_gamma(range));
else
    gamma = mean(sampl_gamma(:,range), 2);
end
R = mean(sampl_R(range));

% Figures for debugging
if ENABLE_DEBUG
    fhl = figure(1);
    set(fhl, 'unit', 'pixel', 'position', [227 299 1241 543]);
    cmap = colormap(lines);
    subplot(3,3,1);
    plot(0:1:iter_max, MHAcceptRatio(1,:), 'o', 'color', cmap(1,:));
    title('Acceptance ratio of sG');
    subplot(3,3,2);
    plot(0:1:iter_max, MHAcceptRatio(2,:), '*', 'color', cmap(2,:));
    title('Acceptance ratio of gamma');
    subplot(3,3,3);
    plot(0:1:iter_max, MHAcceptRatio(3,:), '^', 'color', cmap(3,:));
    title('Acceptance ratio of R');

    subplot(3,3,4);
    plot(0:1:iter_max, MHProbVal(1,:), 'o', 'color', cmap(1,:));
    title('prob of sG');
    subplot(3,3,5);
    plot(0:1:iter_max, MHProbVal(2,:), '*', 'color', cmap(2,:));
    title('prob of gamma');
    subplot(3,3,6);
    plot(0:1:iter_max, MHProbVal(3,:), '^', 'color', cmap(3,:));
    title('prob of R');

    subplot(3,3,7);
    plot(0:1:iter_max, sum(sampl_sG), 'o', 'color', cmap(1,:));
    title('samples of zero-norms of sG');
    subplot(3,3,8);
    if size(sampl_gamma,1) == 1
        plot(0:1:iter_max, sampl_gamma, '*', 'color', cmap(2,:));
    else
        plot(0:1:iter_max, sum(sampl_gamma), '*', 'color', cmap(2,:));
    end
    title('samples of 1-norms of gamma');
    subplot(3,3,9);
    plot(0:1:iter_max, sampl_R, '^', 'color', cmap(3,:));
    title('samples of R');

    fhl2 = figure(2);
    set(fhl2, 'unit', 'pixel', 'position', [288 145 1124 350]);
    for k = 1:size(sampl_sG, 1)
        plot(0:1:iter_max, sampl_sG(k,:)+(k-1)*1.5, '-');
        hold on;
    end
    hold off;
    title('sample path of sG')
end


% ================================================================
% Nested Functions
% ================================================================

function val = evalProb(sG, gamma, R)
% Evaluate marginal probability p(s,gamma,R|t).

% Notes:
% - basic version: diag Gamma, scalar R, so scalar a,c
% - matrix computation: basic, not optimized

    % prior of s: independent Bernoulli distribution
    pval_s = sum(log(pB.^sG + (1-pB).^(1-sG)));

    % prior of gamma and R: Gamma distrition
    if isscalar(gamma)
        gamma = gamma*ones(Ng,1);
    end
    igamma = 1./gamma;
    pval_gamma = sum(log(igamma.^(c-1).*exp(-d*igamma)));
    pval_R = log(R^(a-1)*exp(-b*R));

    % prune Phi
    usePrune = 0;
    if usePrune  % not exact; accelerate computation
        nz_col = find(s);  % nonzero columns in Phi
        Phi_p = Phi(:,nz_col);
        sG_p = sG(nz_col);
        gamma_p = gamma(nz_col);
    else
        Phi_p = Phi;
        sG_p = sG;
        gamma_p = gamma;
    end

    % cond prob of t given s,gamma,R
    s_p = groupLift(sG_p, dimG);
    gamma_p = groupLift(gamma_p, dimG);
    if N < M
        % compute via Sigma_t
        Sigma_t = R*eye(N) + Phi_p*diag(s_p.*gamma_p)*Phi_p';
        pval_t = -.5*log(det(Sigma_t)) - .5*t'*inv(Sigma_t)*t;
    else
        % compute via Sigma_theta
        Sigma_theta_inv = diag(s_p) * Phi_p' * (1/R) * Phi_p * diag(s_p) + ...
            diag(1./gamma_p);
        Sigma_theta = inv(Sigma_theta_inv);
        mu_theta = Sigma_theta*diag(s_p)*Phi_p'/R*t;

        pval_t = -N/2*log(R) -.5*log(prod(gamma_p)) + ...
                 .5*log(det(Sigma_theta)) - ...
                 .5*(t'*t/R - mu_theta'*Sigma_theta_inv*mu_theta);
    end

    val = pval_s + pval_gamma + pval_R + pval_t;
end

end


% ================================================================
% Local Functions
% ================================================================

function vec = findMaxFreq(mat)
% Count the appearance frequency of each element in the given matrix,
% in which each column is an element, and the column index indicates
% the samples.

    % count freq of the 0-norms
    vecSum = sum(mat, 1);
    [freq, bins] = histc(vecSum, unique(vecSum));
    [~, maxFreqIndex] = max(freq);
    [~, indices] = find(bins == maxFreqIndex);
    matRed = mat(:,indices);

    % if all vectors of the given zero-norm are identical, we have found
    % the vector; otherwise, we need to count frequencies further,
    % while now the number of vectors has been significantly reduced.
    [uniqMatRed, im, ium] = unique(matRed', 'rows');
    uniqMatRed = uniqMatRed';
    if size(uniqMatRed,2) == 1
        vec = uniqMatRed;
    else
        % the vectors of the same 0-norm are not identical
        [freq, bins] = histc(ium, unique(ium));
        [~, maxFreqIndex] = max(freq);
        [~, indices] = find(bins == maxFreqIndex);
        vec = matRed(:,indices(1));
    end

end

function vec = groupLift(vecG, dimG)
% This function expands vecG into a full vector according to group
% sizes.

    vec = ones(sum(dimG),1);
    for k = 1:length(dimG)
        if k == 1
            lIndex = 1;
        else
            lIndex = sum(dimG(1:k-1))+1;
        end
        rIndex = sum(dimG(1:k));

        vec(lIndex:rIndex) = vecG(k)*vec(lIndex:rIndex);
    end

end
