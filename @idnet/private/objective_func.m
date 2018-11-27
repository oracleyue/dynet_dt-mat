function optval = objective_func(y, A, w, lambda, grp_midx)
% Calculate the objective function of the optimization problem

ws1norm = 0;
N = size(grp_midx,1);

for k = 1:N
    idxbeg = grp_midx(k,1);
    idxend = grp_midx(k,2);    
    if idxend < idxbeg
        error('Illegal values of the indices of some groups of w_S!')
    end
    
    idxVec = idxbeg:1:idxend;
    sizeGrp = idxend - idxbeg + 1;
    ws1norm = ws1norm + sqrt(sizeGrp) * norm(w(idxVec),2);
end

optval = sumsqr(y - A*w) + lambda * ws1norm;

end