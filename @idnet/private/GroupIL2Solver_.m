function w = GroupIL2Solver_(dy, A, IndexLists, Lambda)
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
% lambda0 = Lambda(3);

% dim information
size_w = size(A,2);

% index maps from w^S1, w^S2 to w
NS1 = C*N; NS2 = N;
for k = 1:NS1
%     if k <= C,  S1Index{k} = (k-1)*blkSize(1)+1 : k*blkSize(1);
    S1Index{k} = C * sum(blkSize(1:ceil(k/C)-1)) * ones(1,blkSize(ceil(k/C))) ...
        + (rem(k-1,C) * blkSize(ceil(k/C))+1 : (rem(k-1,C)+1)*blkSize(ceil(k/C)));
    % S1GroupSize(k) = length(S1Index{k});
    SqrtS1GroupSize(k) = sqrt(length(S1Index{k}));    
end
for k = 1:NS2
    S2Index{k} = C * sum(blkSize(1:k-1)) * ones(1, C*blkSize(k)) ...
        + (1:C*blkSize(k));
    % S2GroupSize(k) = length(S2Index{k});
    SqrtS2GroupSize(k) = sqrt(length(S2Index{k}));
end


% === Optimization
% initialization
epsilon = 1e-4;
w = [];
% iteration
while epsilon >= 1e-8
    % update weights
    if isempty(w)       
        mu = ones(N*C,1);
        nu = ones(N,1);
%         rho = ones(size_w,1);
    else
        mu = (ws1.^2 + epsilon).^(-1/2);
        nu = (ws2.^2 + epsilon).^(-1/2);
%         rho = (w  + epsilon).^(-1);
    end
    % solve reweighted l1
    cvx_begin quiet
        variable w(size_w)
        
        expressions ws1(NS1) ws2(NS2)
        for k = 1:NS1
             ws1(k) = norm(w(S1Index{k}),2);
        end
        for k = 1:NS2
             ws2(k) = norm(w(S2Index{k}),2);
        end
        
        minimize( norm(dy-A*w,2) ...
          + lambda1 * sum_square_pos(mu .* SqrtS1GroupSize' .* ws1) ...
          + lambda2 * sum_square_pos(nu .* SqrtS2GroupSize' .* ws2) )
        %             + lambda0 * norm(diag(rho)*w,1)
    cvx_end

    epsilon = epsilon/10;
end

end