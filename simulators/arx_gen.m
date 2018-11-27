function [A, B, Qb] = arx_gen(num_nodes, spar, order_max, pos_inputs, Ts)
% ARX_GEN generates random ARX models and performs simulations.
% This is a simple version: multiple feedforward, single feedback.
%
% Inputs:
%   - num_nodes: integer, #nodes
%   - spar: (0,1) sparsity density
%   - order_max: 1x2 vector
%       specify the max orders of numerator and denominator
%       polynomials
%   - pos_inputs: i.e. boolean P
%   - Ts: sampling period
% Outputs:
%   - A,B: ARX parametric models: Ay = Bu + e
%   - Qb: network topology

% Note:
% 1. The function "pathfind.m" will cretate a temporary file, named
% "pathTable.csv", which is used to save all paths from node 1 to
% last.
% 2. This is a simple version, which create a digraph with mutiple
% feedforward paths and one feedback path from the last node to the
% first.

% Copyright [2017] <oracleyue>
% Last modified on 18 Jul 2018


% paths
% addpath('./members')    % supportive functions

% parse arguments
if nargin < 4
    pos_inputs = [1; zeros(num_nodes-1,1)];
    Ts = .1;
elseif nargin < 5
    Ts = .1;
end
num_order_max = order_max(1);
den_order_max = order_max(2);

% % data
% num_ts = 1000;          % #points of time series
% % network
% num_nodes = 10;         % #nodes
% spar = .2;              % density of network sparsity
% num_fd = 1;             % number of feedback elements (upper-tri of Q)
% % inputs                % boolean P matrix of DSF
% pos_inputs = [1; zeros(num_nodes-1,1)];
% % element
% den_order_max = 2;      % max order of denominators
% num_order_max = 1;      % max order of numerators
% Ts = .1;                % sampling period

% update density (excluding the feedback path and off-diag paths)
num_ffp = ceil(num_nodes^2*spar) - (num_nodes-1) - 1;
spar_left = num_ffp / (num_nodes^2 - num_nodes) * 2;

% generate toplogy and find all paths
Qb = zeros(num_nodes, num_nodes);
zeroTF = tf(0,1,Ts, 'Variable', 'z');
for ii = 1:num_nodes
    for jj = 1:num_nodes
        if ii > jj
            if ii == jj+1
                Qb(ii,jj) = 1;
            elseif ii > jj+1
                flip_coin = binornd(1,spar_left);
                if flip_coin
                    Qb(ii,jj) = 1;
                else
                    Q(ii,jj) = zeroTF;
                end
            end
        else
            Q(ii,jj) = zeroTF;
        end
    end
end
[out, in] = find(Qb);
if exist('pathTable.csv', 'file')
    delete pathTable.csv
end
pathfind([in'; out'], 1, num_nodes, []);
pathTable = csvread('pathTable.csv');
delete pathTable.csv

% variable to save denominator polynomials
den_cells = cell(num_nodes,1);

% random stable tfs for 1's in topology
[yindex, xindex] = find(Qb');
xidx_old = 0;
for k = 1:length(xindex)
    xidx = xindex(k);
    if xidx_old ~= xidx
        den_poly = stabpoly(den_order_max);
        den_cells{xidx} = den_poly;
    end
    num_order = randi(num_order_max);
    num_poly = stabpoly(num_order);
    tempTF = tf(num_poly, den_poly, Ts, 'Variable', 'z');
    Q(xidx,yindex(k)) = tempTF * (1/hinfnorm(tempTF,1e-4));
    xidx_old = xidx;
end

% -----------
% Note:
% Here is a bug or computational precision in the transfer function
% concatenations. If we compute the total feedforward transfer
% function Tff from node 1 to node LAST, and then compute the H_inf
% norm to design the feeback TF by small gain theorem, it mostly
% failed (getting Inf of Tff) when the size of network
% increases. Instead, we use the norm (triangle) inequality to compute
% the upper bound of H_inf norm of Tff, and then use it to design the
% gain of M.
% -----------

% Method 1:
% use Mason's formula (no loop) to computer the lump-sum TF (Tff) from node
% 1 to node LAST from the above Q (consider all feedforward paths)
%
% for k = 1:size(pathTable,1)
%     path = pathTable(k,:);
%     TPath = tf(1,1,Ts, 'Variable', 'z');
%     ii = 1;
%     while path(ii) < num_nodes
%         TPath = TPath*Q(path(ii+1),path(ii));
%         ii = ii + 1;
%     end
%     if k == 1
%         Tff = TPath;  % tf of all feedforward paths
%     else
%         Tff = Tff + TPath;
%     end
% end

% Method 2:
% instead of computing Tff accurately, we computer the upper bound
% of Hinf norm of Tff, which will be used in small-gain theorem.

% compute Hinf norm of each element of Q
% (no longer needed, since we normalized each element by its Hinf norm)
% gainQ = zeros(size(Q));
% for i = 1:size(Q,1)
%     for j = 1:size(Q,2)
%         if Qb(i,j)
%             gainQ(i,j) = hinfnorm(Q(i,j),1e-4);
%         end
%     end
% end
% (thus, the gainQ is nothing but Qb)
gainQ = Qb;

% compute upper bound of Hinf norm of Tff
gainUpBd = 0;   % upper bound
% use Mason's formula (no loop)
for k = 1:size(pathTable,1)
    path = pathTable(k,:);
    gainPath = 1;
    ii = 1;
    while path(ii) < num_nodes
        gainPath = gainPath * gainQ(path(ii+1),path(ii));
        ii = ii + 1;
    end
    gainUpBd = gainUpBd + gainPath;
end

% adding the stablizable feedback path
num_order = randi(num_order_max);
num_poly = stabpoly(num_order);
den_poly = stabpoly(den_order_max);
den_cells{1} = den_poly;
M = tf(num_poly, den_poly, Ts, 'Variable', 'z');
% turn gain by small-gain theorem
% gainTff = hinfnorm(Tff, 1e-4);
gainTff = gainUpBd;  % use upper bound instead
gainM = hinfnorm(M, 1e-4);
scale = .99/gainTff/gainM;
M = scale*M;
% add feedback path
Q(1,num_nodes) = M;
Qb(1,num_nodes) = 1;


% generate P matrix ([1; 0; 0; ... 0])
Pb = pos_inputs;
for k = 1:size(Pb,1)
    den_poly = den_cells{k};
    for l = 1:size(Pb,2)
        if Pb(k,l)
            num_order = randi(num_order_max);
            num_poly = stabpoly(num_order);
            P(k,l) = tf(num_poly, den_poly, Ts, 'Variable', 'z');
        else
            P(k,l) = zeroTF;
        end
    end
end


% convert to ARX models

% % model format 1: (/Yue, IFAC 2017/)
% A = den_cells;
% [By, ~] = tfdata(Q);
% [Bu, ~] = tfdata(P);

% model format 2:
[A, ~] = tfdata(Q);
[B, ~] = tfdata(P);
for k = 1:num_nodes
    A{k,k} = den_cells{k};
end
