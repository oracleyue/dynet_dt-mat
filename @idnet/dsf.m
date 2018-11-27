function dynetObj = dsf(idnetObj)
% Convert idnet-Object to dynet-Object
%
% Input:
%   idnetObj: idnet object
%        - A : px1 cell, A{i} = [1, a_{i1}, ..., a_{i n_i^a}], diagonal
%        - By: pxp cell, By{i,j} = [b^y_{ij1}, ..., b^y_{ij n_ij^by}]
%        - Bu: pxp cell, Bu{i,j} = [b^u_{ij1}, ..., b^u_{ij n_ij^bu}]
%        - C : px1 cell, C{i} = [1, c_{i1}, ..., c_{i n_i^c}], diagonal
%        - model: A(q)y(t) = By(q)y(t) + Bu(q)u(t) + C(q)e(t)
%
% Output:
%   dynetObj: dynet object, structure of {Q, P, H, Ts}
%        - Q : pxp tf
%        - P : pxm tf
%        - H : px1 tf, diagnoal
%        - Ts: double
%        - model: y(t) = Q(q)y(t) + P(q)u(t) + H(q)e(t)
%

    numsQ = addZeroToB(idnetObj.By);
    numsP = addZeroToB(idnetObj.Bu);
    densQ = cellConcatenate(idnetObj.A, size(idnetObj.By,2));    
    densP = cellConcatenate(idnetObj.A, size(idnetObj.Bu,2));
    numsH = idnetObj.C;
    densH = idnetObj.A;
    Ts = idnetObj.Ts;
    
    dynetObj.Q = idtf(numsQ, densQ, Ts);
    dynetObj.P = idtf(numsP, densP, Ts);
    if isempty(numsH)
        dynetObj.H = idtf(numsH, densH, Ts);        
    end
    dynetObj.Ts = idnetObj.Ts;
    dynetObj.Variable = 'z^-1';
end


% ====================================
% Local functions
% ====================================
function AA = cellConcatenate(A, p)    
    AA = {};
    for i = 1:p
        AA = [AA A'];
    end
end

function B = addZeroToB(B)
    [n,m] = size(B);
    for i = 1:n
        for j = 1:m
            B{i,j} = [0 B{i,j}];
        end
    end    
end