function w = GroupL1Solver_(by, A, IndexLists, Lambda, config)
% Solve "dy = A wi + lambda1 ||wi||_0^S1 + lambda2 ||wi||_0^S2"
%
% Inputs:
%     IndexLists = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
%     Lambda = [lambda1, lambda2]

% settings for the solver
glb_optimizer = config.name;
MAX_ITER = config.max_iter;
ABSTOL = config.abstol;

% unpack input variables    
N = IndexLists{1}(1);
C = IndexLists{1}(2);
blksize = IndexLists{2};
lambda = Lambda;

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
        size_w = size(A,2);

        cvx_begin quiet
            variable w(size_w) 

            expressions ws(N)
            for k = 1:N
                idxVec = grp_midx(k,1):1:grp_midx(k,2);
                ws(k) = sqrt(grp_size(k)) * norm(w(idxVec),2);
            end

            minimize( norm(by-A*w,2) + lambda * sum(ws) )    
        cvx_end       
        % == END: CVX
    
    case 'prox'
        % define precisions        
%         global MAX_ITER ABSTOL
        if isempty(MAX_ITER) || isempty(ABSTOL)
            MAX_ITER = 1000;
            ABSTOL = 1e-6;
        end
        % == Optimization via Proximal Gradient ==  
        [w, ~] = proxgrad(by, A, lambda, grp_midx, [MAX_ITER ABSTOL]);
              
end
