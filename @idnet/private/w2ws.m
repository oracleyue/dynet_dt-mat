function ws = w2ws(w, grp_midx)
% To obtain ws from w

N = size(grp_midx,1);
ws = zeros(N,1);

for k = 1:N   
    idxbeg = grp_midx(k,1);
    idxend = grp_midx(k,2);    
    if idxend < idxbeg
        error('Illegal values of the indices of some groups of w_S!')
    end
    
    idxVec = idxbeg:1:idxend;
    sizeGrp = idxend - idxbeg + 1;
    ws(k) = norm(w(idxVec),2);
end

end