function [w, pred_err] = NGroupIL1Solver_(by, A, IndexList, LambdaList, config)
% Solve "dy = A wi + lambda_WE ||wi^E||_0 + lambda_WS ||wi^S||_0"
%
% Inputs:
%     IndexList = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
%     LambdaList = [lambda_WE, lambda_WS, inf]
%
% Outputs:
%     pred_err: it could be opt_val, rmse, cod_k
%
% Outputs (depricated):
%     cod_k, rmse


% settings for the solver
glb_optimizer = config.name;
MAX_ITER = config.max_iter;
ABSTOL = config.abstol;

% unpack input variables
N = IndexList{1}(1);
C = IndexList{1}(2);
blksize = IndexList{2}';
[Ny, Nw] = size(A);
% lambda = LambdaList(2);
lambda_wE = LambdaList(1);
lambda_wS = LambdaList(2);

NE = C*N; NS = N;
% index map of small group: w^E
grp_size_wE = kron(blksize, ones(C,1));
grp_midx_wE = zeros(NE, 2);
for k = 1:NE
    % hint: if k <= C,  tmpidx = (k-1)*blksize(1)+1 : k*blksize(1);
    tmpidx = C * sum(blksize(1:ceil(k/C)-1)) * ones(1,blksize(ceil(k/C))) ...
        + (rem(k-1,C) * blksize(ceil(k/C))+1 : (rem(k-1,C)+1)*blksize(ceil(k/C)));
    grp_midx_wE(k,:) = [tmpidx(1) tmpidx(end)];
end

% index map of big group: w^S
grp_size_wS = C*blksize;
grp_midx_wS = zeros(NS, 2);
for k = 1:NS
    tmpidx = C * sum(blksize(1:k-1)) * ones(1, C*blksize(k)) ...
             + (1:C*blksize(k));
    grp_midx_wS(k,:) = [tmpidx(1) tmpidx(end)];
end

% choose the solver for optimizaiton
switch glb_optimizer
    case 'cvx'
        % == Optimization via CVX ==
        rhoE = sqrt(grp_size_wE);
        rhoS = sqrt(grp_size_wS);
        epsilon = 1e-4;   % initialize \epsilon
        w = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                mu = ones(NE,1);
                nu = ones(NS,1);
            else
                mu = (wE + epsilon).^(-1);
                nu = (wS + epsilon).^(-1);
            end
            % solve reweighted l1
            cvx_begin quiet
                variable w(Nw)
                expressions wE(NE) wS(NS)
                    for k = 1:NE
                        idxVec = grp_midx_wE(k,1):1:grp_midx_wE(k,2);
                        wE(k) = sqrt(grp_size_wE(k)) * norm(w(idxVec),2);
                    end
                    for k = 1:NS
                        idxVec = grp_midx_wS(k,1):1:grp_midx_wS(k,2);
                        wS(k) = sqrt(grp_size_wS(k)) * norm(w(idxVec),2);
                    end
                minimize( norm(by-A*w,2) + lambda_wE*sum(mu .* rhoE .* wE) ...
                                         + lambda_wS*sum(nu .* rhoS .* wS) )
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
        w = []; wS = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (wS + epsilon).^(-1);
            end
            % solve reweighted-l1 optimization
            [w, ~] = proxgrad_wt(by, A, lambda, grp_midx_wS, nu, [MAX_ITER ABSTOL]);
            wS = w2ws(w, grp_midx_wS);
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
        w = []; wS = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (wS + epsilon).^(-1);
            end
            % solve reweighted-l1 optimization
            [w, ~] = proxgradacc_wt(by, A, lambda, grp_midx_wS, nu, [MAX_ITER ABSTOL]);
            wS = w2ws(w, grp_midx_wS);
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
        w = []; wS = [];
        % iteration
        while epsilon >= 1e-8
            % update weights
            if isempty(w)
                nu = ones(N,1);
            else
                nu = (wS + epsilon).^(-1);
            end
            % solve reweighted-l1 optimization
            [w, ~] = proxadmm_wt(by, A, lambda, grp_midx_wS, nu, [MAX_ITER ABSTOL]);
            wS = w2ws(w, grp_midx_wS);
            % update \epsilon
            epsilon = epsilon/10;
        end
        % == END: ADMM
end
% cod_k = cod_calc(by, A, w, 1);
% rmse = rmse_calc(by, A, w);

% pred_err = cod_calc(by, A, w, 1);
% pred_err = rmse_calc(by, A, w);
pred_err = objective_func_nl(by, A, w, lambda_wE, lambda_wS,...
                             grp_midx_wE, grp_midx_wS);

end