function optval = objective_func_nl(y, A, w, lambda_wE, lambda_wS, ...
                                    grp_midx_wE, grp_midx_wS)
% Calculate the objective function of the optimization problem
% version for nonlinear dynamic networks

wsEnorm = 0;
wsSnorm = 0;
NE = size(grp_midx_wE,1);
NS = size(grp_midx_wS,1);

for k = 1:NE
    idxbeg = grp_midx_wE(k,1);
    idxend = grp_midx_wE(k,2);
    if idxend < idxbeg
        error('Illegal values of the indices of some groups of w_S!')
    end

    idxVec = idxbeg:1:idxend;
    sizeGrp = idxend - idxbeg + 1;
    wsEnorm = wsEnorm + sqrt(sizeGrp) * norm(w(idxVec),2);
end

for k = 1:NS
    idxbeg = grp_midx_wS(k,1);
    idxend = grp_midx_wS(k,2);
    if idxend < idxbeg
        error('Illegal values of the indices of some groups of w_S!')
    end

    idxVec = idxbeg:1:idxend;
    sizeGrp = idxend - idxbeg + 1;
    wsSnorm = wsSnorm + sqrt(sizeGrp) * norm(w(idxVec),2);
end

optval = sumsqr(y - A*w) + lambda_wE * wsEnorm + lambda_wS * wsSnorm;

end