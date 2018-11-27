function [w, pred_err] = GroupIL1Solver_(by, A, IndexLists, optimizerOpt)
% Solve "dy = A wi + lambda1 ||wi||_0^S1 + lambda2 ||wi||_0^S2"
%
% Inputs:
%     IndexLists = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
%     Lambda = [lambda1, lambda2]
% Outputs:
%     pred_err: it could be opt_val, rmse, cod_k
% Outputs (depricated):
%     cod_k, rmse

% settings for the solver

glb_optimizer = optimizerOpt.method;
MAX_ITER = optimizerOpt.MAXITER;
ABSTOL = optimizerOpt.ABSTOL;

% unpack input variables
N = IndexLists{1}(1);
C = IndexLists{1}(2);
blksize = IndexLists{2};
[Ny, Nw] = size(A);
lambda = optimizerOpt.lambda;

% grouping
grp_size = C*blksize;
grp_midx = zeros(N, 2);
for k = 1:N
    tmpidx = C * sum(blksize(1:k-1)) * ones(1, C*blksize(k)) ...
             + (1:C*blksize(k));
    grp_midx(k,:) = [tmpidx(1) tmpidx(end)];
end

% choose the solver for optimizaiton
switch glb_optimizer
    case 'cvx'
        % == Optimization via CVX ==
        rhoS = sqrt(grp_size);
        epsilon = 1e-4;   % initialize \epsilon
        w = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (ws + epsilon).^(-1);
            end
            % solve reweighted l1
            cvx_begin quiet
                variable w(Nw)
                expressions ws(N)
                    for k = 1:N
                        idxVec = grp_midx(k,1):1:grp_midx(k,2);
                        ws(k) = sqrt(grp_size(k)) * norm(w(idxVec),2);
                    end
                minimize( norm(by-A*w,2) + lambda*sum(nu .* rhoS' .* ws) )
            cvx_end
            % update \epsilon
            epsilon = epsilon/10;
        end
        % == END: CVX


    case 'prox'
        % == Optimization via Proximal Gradient Method ===
%         prox_timer = tic;
        % define precisions
%         global MAX_ITER ABSTOL
        if isempty(MAX_ITER) || isempty(ABSTOL)
            MAX_ITER = 1e6;
            ABSTOL = 1e-12;
        end
        % initialize \epsilon
        epsilon = 1e-4;
        w = []; ws = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (ws + epsilon).^(-1);
            end
            % solve reweighted-l1 optimization
            [w, ~] = proxgrad_wt(by, A, lambda, grp_midx, nu, [MAX_ITER ABSTOL]);
            ws = w2ws(w, grp_midx);
            % update \epsilon
            epsilon = epsilon/10;
        end
%         toc(prox_timer)
        % == END: Proximal Gradient Method


    case 'proxacc'
        % == Optimization via Accelerated Proximal Gradient Method ===
        if isempty(MAX_ITER) || isempty(ABSTOL)
            MAX_ITER = 1e6;
            ABSTOL = 1e-12;
        end
        % initialize \epsilon
        epsilon = 1e-4;
        w = []; ws = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (ws + epsilon).^(-1);
            end
            % solve reweighted-l1 optimization
            [w, ~] = proxgradacc_wt(by, A, lambda, grp_midx, nu, [MAX_ITER ABSTOL]);
            ws = w2ws(w, grp_midx);
            % update \epsilon
            epsilon = epsilon/10;
        end
        % == END: Proximal Gradient Method


    case 'admm'
        % == Optimization via ADMM ===
        if isempty(MAX_ITER) || isempty(ABSTOL)
            MAX_ITER = 1e6;
            ABSTOL = 1e-12;
        end
        % initialize \epsilon
        epsilon = 1e-4;
        w = []; ws = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (ws + epsilon).^(-1);
            end
            % solve reweighted-l1 optimization
            [w, ~] = proxadmm_wt(by, A, lambda, grp_midx, nu, [MAX_ITER ABSTOL]);
            ws = w2ws(w, grp_midx);
            % update \epsilon
            epsilon = epsilon/10;
        end
        % == END: ADMM
end
% cod_k = cod_calc(by, A, w, 1);
% rmse = rmse_calc(by, A, w);

% pred_err = cod_calc(by, A, w, 1);
% pred_err = rmse_calc(by, A, w);
pred_err = objective_func(by, A, w, lambda, grp_midx);

end