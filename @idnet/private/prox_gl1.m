function prox = prox_gl1(v, size, gam_lam)
% Calculate the proximal opeartor of group l1-norm
% gam_lam = gamma * lambda
    
radius = gam_lam * sqrt(size);
v2norm = norm(v,2);
prox = (1 - radius/v2norm) * v * (v2norm >= radius);

end