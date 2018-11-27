function w = IL2Solver_(dy, A, IndexLists, Lambda)
% Solve "dy = A wi + lambda1 ||wi||_0^S1 + lambda2 ||wi||_0^S2"
%
% Inputs:
%     IndexLists = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
%     Lambda = [lambda1, lambda2]


% unpack input variables    
N = IndexLists{1}(1);
C = IndexLists{1}(2);
blkSize = IndexLists{2};
lambda1 = Lambda(1);
lambda2 = Lambda(2);
lambda0 = Lambda(3);

% dim information
size_w = size(A,2);


% === Optimization
% initialization
epsilon = 1e-4;
w = [];
% iteration
while epsilon >= 1e-8
    % update weights
    if isempty(w)       
        rho = ones(size_w,1);
    else
%         rho = (w.^2 + epsilon).^(-1/2);
        rho = (w.^2 + epsilon).^(-1/2);
    end
    % solve reweighted l2
    cvx_begin quiet
        variable w(size_w)
        expression ws
        for k = 1:size_w
            ws(k) = rho(k) * norm(w(k),2);
        end
        minimize( norm(dy-A*w,2) ...
              + lambda0 * sum_square_pos(ws) )
%             + lambda0 * sum_square(diag(rho)*w) )
    cvx_end

    epsilon = epsilon/10;
end

end