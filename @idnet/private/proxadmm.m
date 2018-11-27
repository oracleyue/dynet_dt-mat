function [w, solver_info] = proxadmm(by, A, lambda, grp_midx, precision)
% Accelearated proxgrad, equipped with Line Search
%
% Input:
%   - by, A  : (M x 1) vector, (M x size_w) matrix
%              by = A*w + noise
%   - lambda : double
%              regularization parameters
%   - nu : (N x 1) weight vector
%          weights for L1-reweighted method    
%   - precision: [MAX_ITER ABSTOL]
%        
% 
% Output:
%   - w      : (size_w x 1) vector
%   - solver_info : structure
%         the status of the optimization solution     
%

% precisions
MAX_ITER = precision(1);
ABSTOL = precision(2);

% dataset parameters
N = size(grp_midx,1);
size_w = size(A,2);
grp_size = grp_midx(:,2) - grp_midx(:,1) + ones(N,1);

% dimension checking
if length(nu) ~= N
    error('The weight vector \nu fails to match the group structure.')
end


% == Optimization via Accelerated Proximal Gradient Method ===

proxacc_timer = tic;    % /2/tic

% precompute
AtA = A'*A;
Aty = A'*by;
L = svds(A,1);
gamma = 1/L;

% declare variables and temporary functions
size_w = size(A,2);
w = zeros(size_w,1);
admm_status = 'Unsolved';
z = w; u = w;
I = eye(size_w);
% ... for line search 
wtry = w;
f = @(w) sumsqr(A*w-by);
beta = .5;

% initialize variables
% admm_optval = zeros(MAX_ITER, 1);
admm_optval = zeros(2,1);

% iteration without /line search/
for l = 1:MAX_ITER
    while 1   % Line Search        
        wtry = (I + gamma*AtA)\(gamma*Aty + z - u);        
        for k = 1:N
            idxVec = grp_midx(k,1):1:grp_midx(k,2);
            vk = wtry(idxVec) + u(idxVec);
            proxgk = prox_gl1(vk, grp_size(k), gamma*lambda); 
            z(idxVec) = proxgk;
        end
        u = u + wtry - z;
        if f(wtry) <= f(w) + 2*u'*(wtry-w) + (1/2/gamma)*sumsquare(wtry-w)
            break;
        end
        gamma = beta * gamma;
    end
    w = wtry;
        
    % admm_optval(l) = objective_func(by, A, w, lambda, grp_midx);    
    admm_optval(1) = admm_optval(2);
    admm_optval(2) = objective_func(by, A, w, lambda, grp_midx);
    
    if l > 1 && abs(admm_optval(2) - admm_optval(1)) < ABSTOL
    % if l > 1 && abs(admm_optval(l) - admm_optval(l-1)) < ABSTOL
%             /abs(admm_optval(l)) < RELTOL
        admm_status = 'Solved';
        break;
    end
end


% Prepare status results
solver_info.status = admm_status;
solver_info.p = admm_optval(end);
solver_info.cputime = toc(proxacc_timer);
solver_info.iterates = l;

% admm_cputime = toc   % /2/toc
% admm_w = w;
% admm_p = admm_optval(end)
% admm_iter = l;

end