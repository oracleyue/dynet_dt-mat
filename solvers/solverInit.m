function solverOptions = solverInit(solver, lambda, learnLambda)
% SOLVERINIT initialize the options for solvers used in network
% inference.

% OPTIONAL INPUT:
%   solver: string, "gsbl", "gil1", "gl1"
%   lambda: double OR [double, double, double] (see "arx.m")
%           regularization parameter for GLASSO or initial value of
%           noise variance for GSBL or GSMC
%   learnLambda: 1 or 0

% OUTPUT:
%   solverOptions: struct
%      - solver: string, "gl1", "gil1", "gsbl", "l1", "il1"; by default "gil1"
%      - lambda: double OR [double, double, double] (see "arx.m")
%      - optimizerOpt: struct, to control numerical optimization
%           - method: string, "cvx", "prox", "proxacc", "admm" (unsed for "gsbl")
%           - LearnLambda: 1 or 0  (only for "gsbl)
%           - MAXITER: maximal iterations
%           - ABSTOL: absolute tolerance
%           - RELTOL: relative tolerance  (unsed for "gsbl")
%           - PruneGamma: threshold to prune gamma (only for "gsbl)

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 22 Aug 2017


% arugment parsing
switch nargin
  case 0
    solverOptions.solver = 'gsbl';
  otherwise
    solverOptions.solver = solver;
end

% setting default values
switch solverOptions.solver
  case 'gsbl'
    switch nargin
      case 1
        solverOptions.noiseVar0 = 1e-2;
        solverOptions.LearnLambda = 1;
      case 2
        solverOptions.noiseVar0 = lambda;
        solverOptions.LearnLambda = 1;
      case 3
        solverOptions.LearnLambda = learnLambda;
    end
    solverOptions.method      = 'EM';
    solverOptions.MAXITER     = 2000;
    solverOptions.ABSTOL      = 1e-8;
    solverOptions.PruneGamma  = 1e-3;

  case 'gil1'
    switch nargin
      case 1
        solverOptions.lambda = 1e-2;
      otherwise
        solverOptions.lambda = lambda;
    end
    solverOptions.method  = 'cvx';
    solverOptions.MAXITER = 1e4;    % used for "prox"
    solverOptions.ABSTOL  = 1e-12;  % used for "prox"
    solverOptions.RELTOL  = 1e-8;   % used for "prox"

  case 'gspmc'
    hypara.a = 1e-4;
    hypara.b = 1e-4;
    hypara.c = 1e-4;
    hypara.d = 1e-4;
    hypara.pB = .5;
    switch nargin
      case 1
        hypara.R0 = 1e-2;
      otherwise
        hypara.R0 = lambda;
    end
    hypara.gamma0 = 1;

    mcpara.iter_max = 4e3;  % max number of Gibbs iterations
    mcpara.var_gamma = .01;
    mcpara.var_R = .01;
    mcpara.prunePhi = 0;
    mcpara.MH_scheme = 1;  % 1 or 2; schemes of proposal dist of s
    mcpara.Pred_scheme = 'mean';
    mcpara.ENABLE_DEBUG = 0;

    solverOptions.method = 'Gibbs M-H';
    solverOptions.hypara = hypara;
    solverOptions.mcpara = mcpara;
    solverOptions.confid = .95;
    % confidence, i.e., val := 0 if val < .05; val := 1 if val > .95
end
