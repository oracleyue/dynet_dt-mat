function [w, solver_info] = proxgradacc(by, A, lambda, grp_midx, precision)
% Accelearated proxgrad, equipped with Line Search
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

% precisions
MAX_ITER = precision(1);
ABSTOL = precision(2);

% dataset parameters
N = size(grp_midx,1);
size_w = size(A,2);
grp_size = grp_midx(:,2) - grp_midx(:,1) + ones(N,1);


% == Optimization via Accelerated Proximal Gradient Method ===

proxacc_timer = tic;    % /2/tic

% precompute
AtA = A'*A;
Aty = A'*by;
L = svds(A,1);
gamma = 1/L;

% declare variables and temporary functions
w = zeros(size(A,2),1);
prox_status = 'Unsolved';
wprev = w;  x = w;
% ... for line search 
z = w;
f = @(w) sumsqr(A*w-by);
beta = .5;

% initialize variables
% prox_optval = [];
prox_optval = zeros(2,1);

% iteration without /line search/
for l = 1:MAX_ITER
    while 1   % Line Search
        mu = l/(l+3);
        x = w + mu*(w - wprev);
        u = AtA*x - Aty;        
        for k = 1:N
            idxVec = grp_midx(k,1):1:grp_midx(k,2);
            vk = x(idxVec) - gamma*u(idxVec);
            proxk = prox_gl1(vk, grp_size(k), gamma*lambda);                
            z(idxVec) = proxk;
        end                
        if f(z) <= f(w) + 2*u'*(z-w) + (1/2/gamma)*sumsquare(z-w)
            break;
        end
        gamma = beta * gamma;
    end
    wprev = w;
    w = z;
        
    % prox_optval(l) = objective_func(by, A, w, lambda, grp_midx);    
    prox_optval(1) = prox_optval(2);
    prox_optval(2) = objective_func(by, A, w, lambda, grp_midx);
    
    if l > 1 && abs(prox_optval(2) - prox_optval(1)) < ABSTOL    
    % if l > 1 && abs(prox_optval(l) - prox_optval(l-1)) < ABSTOL
%             /abs(prox_optval(l)) < RELTOL
        prox_status = 'Solved';
        break;
    end
end


% Prepare status results
solver_info.status = prox_status;
solver_info.p = prox_optval(end);
solver_info.cputime = toc(proxacc_timer);
solver_info.iterates = l;

% prox_cputime = toc   % /2/toc
% prox_w = w;
% prox_p = prox_optval(end)
% prox_iter = l;

end