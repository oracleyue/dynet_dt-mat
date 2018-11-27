function [w sigma2] = ElementSBLSolver_(dy, A, indexLists, lambda, optimizerOpt)
% Solve "dy = A wi + lambda1 ||wi||_0^S1 + lambda2 ||wi||_0^S2"
%
% Inputs:
%     indexLists = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
%     lambda
% Ouputs:
%     w: estimation of parameters
%     sigma2: estimation of unknown noise variance (\Sigma = sigma2 I)

% Copyright [2017] <oracleyue>
% Last modified on 23 Jul 2018


%
% unpack arguments
%
N = indexLists{1}(1);
C = indexLists{1}(2);
blksize = indexLists{2};    % row vector
[Ny, Nw] = size(A);


%
% build up group index table
%
grp_size = C*blksize;
grp_midx = zeros(N, 2);
grp_midx = [cumsum([1 grp_size(1:end-1)]); cumsum(grp_size)]';


%
% List of solvers for elementary sparsity
%

% - use solvers from Wipf
% [w,gamma_est,gamma_used,count] = MSBL(A, dy, lambda, 1);

% - use solvers modified by Yue
[w, sigma2, gamma_est, SBL_SolverStatus] = ...
    classicSBL(A, dy, grp_midx, lambda, optimizerOpt);


end
