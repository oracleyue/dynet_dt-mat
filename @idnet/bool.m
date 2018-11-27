function boolnet = bool(net, dynet_type, zero)
% To generate boolean structure from idnet object.
%
% Ouptut:
%   boolnet: a structure of two entries "Q" and "P"
%       - boolnet.Q: a pxp matrix of 1-0
%       - boolnet.P: a pxm matrix of 1-0
%
%
% Info of idnet object:
% - Linear Dynamic Networks:
%   A : px1 cell, A{i} = [1, a_{i1}, ..., a_{i n_i^a}], diagonal
%   By: pxp cell, By{i,j} = [b^y_{ij1}, ..., b^y_{ij n_ij^by}]
%   Bu: pxp cell, Bu{i,j} = [b^u_{ij1}, ..., b^u_{ij n_ij^bu}]
%   C : px1 cell, C{i} = [1, c_{i1}, ..., c_{i n_i^c}], diagonal
% - Nonlinear Boolean Dynamic Networks:
%   A : px1xC cell, A{:,:,c} is the A for the c-th experiment
%   By: pxpxC cell
%   Bu: pxpxC cell
%   C : px1xC cell
%
%   - model: A(q)y(t) = By(q)y(t) + Bu(q)u(t) + C(q)e(t)

if nargin == 1
    zero = net.ThresholdZero;
    dynet_type = 'L';
elseif nargin == 2
    zero = net.ThresholdZero;
end

p = size(net.By, 1);
m = size(net.Bu, 2);
if iscell(net.By{1,2})
    C = length(net.By{1,2});
    type = 'multiple';  % multiple experiments/datasets
else
    type = 'single';    % single experiment
end

switch dynet_type
  % ------------------------------------------------
  case {'L', 'l'}   % Linear Dynamic Networks

    switch type
      case 'multiple'
        for c = 1:C
            tempBoolQ = zeros(p,p);
            tempBoolP = zeros(p,m);
            for i = 1:p
                for j = 1:p
                    if norm(net.By{i,j}{c}, 2) <= zero
                        tempBoolQ(i,j) = 0;
                    else
                        tempBoolQ(i,j) = 1;
                    end
                end
                for k = 1:m
                    if norm(net.Bu{i,k}{c}, 2) <= zero
                        tempBoolP(i,k) = 0;
                    else
                        tempBoolP(i,k) = 1;
                    end
                end
            end

            if c == 1
                boolQ = tempBoolQ;
                boolP = tempBoolP;
            else
                boolQ = tempBoolQ & boolQ;
                boolP = tempBoolP & boolP;
            end
        end

      case 'single'
        tempBoolQ = zeros(p,p);
        tempBoolP = zeros(p,m);
        for i = 1:p
            for j = 1:p
                if norm(net.By{i,j}, 2) <= zero
                    tempBoolQ(i,j) = 0;
                else
                    tempBoolQ(i,j) = 1;
                end
            end
            for k = 1:m
                if norm(net.Bu{i,k}, 2) <= zero
                    tempBoolP(i,k) = 0;
                else
                    tempBoolP(i,k) = 1;
                end
            end
        end
    end

    boolnet.Q = tempBoolQ;
    if m % having inputs
        boolnet.P = tempBoolP;
    end

  % ------------------------------------------------
  case {'NL', 'nl'}  % Boolean Nonlinear Dynamic Networks

    for c = 1:C
        tempBoolQ = zeros(p,p);
        tempBoolP = zeros(p,m);
        for i = 1:p
            for j = 1:p
                if norm(net.By{i,j}{c}, 2) <= zero
                    tempBoolQ(i,j) = 0;
                else
                    tempBoolQ(i,j) = 1;
                end
            end
            for k = 1:m
                if norm(net.Bu{i,k}{c}, 2) <= zero
                    tempBoolP(i,k) = 0;
                else
                    tempBoolP(i,k) = 1;
                end
            end
        end

        if c == 1
            boolQ = tempBoolQ;
            boolP = tempBoolP;
        else
            boolQ = tempBoolQ | boolQ;
            boolP = tempBoolP | boolP;
        end
    end

    boolnet.Q = boolQ;
    if m % having inputs
        boolnet.P = boolP;
    end

  % ------------------------------------------------
  case {'NLL', 'nll'}  % Boolean Nonlinear Dynamic Networks; Save All Bool Subnetworks

    for c = 1:C
        tempBoolQ = zeros(p,p);
        tempBoolP = zeros(p,m);
        for i = 1:p
            for j = 1:p
                if norm(net.By{i,j}{c}, 2) <= zero
                    tempBoolQ(i,j) = 0;
                else
                    tempBoolQ(i,j) = 1;
                end
            end
            for k = 1:m
                if norm(net.Bu{i,k}{c}, 2) <= zero
                    tempBoolP(i,k) = 0;
                else
                    tempBoolP(i,k) = 1;
                end
            end
        end

        if c == 1
            boolQ = tempBoolQ;
            boolP = tempBoolP;
        else
            boolQ = tempBoolQ | boolQ;
            boolP = tempBoolP | boolP;
        end

        boolQArray{c} = tempBoolQ;
        boolPArray{c} = tempBoolP;
    end

    boolnet.Q = boolQ;
    boolnet.QArray = boolQArray;
    if m % having inputs
        boolnet.P = boolP;
        boolnet.PArray = boolPArray;
    end

end