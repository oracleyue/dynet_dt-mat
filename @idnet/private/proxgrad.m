function [w, solver_info] = proxgrad(by, A, lambda, grp_midx, precision)
% Proximal Gradient Method with Line Search to solve grouped L1
%
% Input:
%   - by, A  : (M x 1) vector, (M x size_w) matrix
%              by = A*w + noise
%   - lambda : double
%              regularization parameters
%   - precision: [MAX_ITER ABSTOL]
%        
% 
% Output:
%   - w      : (size_w x 1) vector
%   - solver_info : structure
%         the status of the optimization solution     
%


% % declare global variables
% global MAX_ITER ABSTOL   
% if isempty(MAX_ITER) || isempty(ABSTOL)
%     error('The parameters "MAX_ITER, ABSTOL" for precision are not defined!')
% end
MAX_ITER = precision(1);
ABSTOL = precision(2);

% dataset parameters
N = size(grp_midx,1);
% N = length(blksize);
size_w = size(A,2);
% C = size_w/N  % Wrong! has to be passed via func parameters
grp_size = grp_midx(:,2) - grp_midx(:,1) + ones(N,1);

% == Optimization via Proximal Gradient Method ===

% % grouping 
% grp_size = C*blksize;
% grp_midx = zeros(N, 2);
% for k = 1:N
%     tmpidx = C * sum(blksize(1:k-1)) * ones(1, C*blksize(k)) ...
%              + (1:C*blksize(k));   
%     grp_midx(k,:) = [tmpidx(1) tmpidx(end)];
% end

prox_timer = tic;

% precompute
AtA = A'*A;
Aty = A'*by;
L = svds(A,1);
gamma = 1/L;

% precisions
% MAX_ITER = 1000;
% ABSTOL = 1e-6;
% RELTOL = 1e-4;

% declare variables and temporary functions
w = zeros(size(A,2),1);
prox_status = 'Unsolved';
% ... for line search 
z = w;
f = @(w) sumsqr(A*w-by);
beta = .5;

% initialize variables
prox_optval = [];

% iteration without /line search/
for l = 1:MAX_ITER
    while 1   % Line Search
        u = AtA*w - Aty;  % gradient
        for k = 1:N
            idxVec = grp_midx(k,1):1:grp_midx(k,2);
            uk = u(idxVec);
            vk = w(idxVec) - gamma*uk;
            proxk = prox_gl1(vk, grp_size(k), gamma*lambda);                
            z(idxVec) = proxk;
        end                
        if f(z) <= f(w) + 2*u'*(z-w) + (1/2/gamma)*sumsquare(z-w)
            break;
        end
        gamma = beta * gamma;
    end
    w = z;    
    
    prox_optval(l) = objective_func(by, A, w, lambda, grp_midx);
    
    if l > 1 && abs(prox_optval(l) - prox_optval(l-1)) < ABSTOL
%             /abs(prox_optval(l)) < RELTOL
        prox_status = 'Solved';
        break;
    end
end

% Prepare status results
solver_info.status = prox_status;
solver_info.p = prox_optval(end);
solver_info.cputime = toc(prox_timer);
solver_info.iterates = l;

% prox_cputime = toc(prox_timer);
% prox_w = w;
% prox_p = prox_optval(end);
% prox_iter = l;



end