function net = ErrorEdgeReduction(Ggt, net, s,t, val)

N = length(s);

for i = 1:N
    sp = s(i); tp = t(i);
    path = shortestpath(Ggt, sp, tp);
    
    if length(path) ~= 0 && length(path) <= 2
       net(sp, tp) = 0;
    end    
end

end
