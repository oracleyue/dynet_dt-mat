function [w sigma2] = GroupSBLSolver_(dy, A, indexLists, optimizerOpt)
% Solve "dy = A wi + lambda1 ||wi||_0^S1 + lambda2 ||wi||_0^S2"
%
% Inputs:
%     indexLists = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
%     lambda
% Ouputs:
%     w: estimation of parameters
%     sigma2: estimation of unknown noise variance (\Sigma = sigma2 I)

% Copyright [2017] <oracleyue>
% Last modified on 22 Aug 2017


%
% unpack arguments
%
N = indexLists{1}(1);
C = indexLists{1}(2);
blksize = indexLists{2};    % row vector
[Ny, Nw] = size(A);

lambda = optimizerOpt.noiseVar0;


%
% build up group index table
%
grp_size = C*blksize;
grp_midx = zeros(N, 2);
% for k = 1:N
%     tmpidx = C * sum(blksize(1:k-1)) * ones(1, C*blksize(k)) ...
%              + (1:C*blksize(k));
%     grp_midx(k,:) = [tmpidx(1) tmpidx(end)];
% end

% alternative to build group index table
grp_midx = [cumsum([1 grp_size(1:end-1)]); cumsum(grp_size)]';


%
% (test) List of solvers for elementary sparsity
%

% - use solvers from Wipf
% [w,gamma_est,gamma_used,count] = MSBL(A, dy, lambda, 1);

% - use solvers from Tipping
% OPTIONS		= SB2_UserOptions('iterations', 500, ...
%                               'fixedNoise', true, ...
% 							  'diagnosticLevel', 0);
% SETTINGS	= SB2_ParameterSettings('NoiseStd', lambda);

% [parameter, hyperParameter, diagnostic] = ...
%     SparseBayes('Gaussian', A, dy, OPTIONS, SETTINGS);

% w = zeros(size_w,1);
% w(parameter.Relevant)	= parameter.Value;



%
% List of solvers for group sparsity
%
[w, sigma2, gamma_est, SBL_SolverStatus] = ...
    groupSBL(A, dy, grp_midx, lambda, optimizerOpt);


end
