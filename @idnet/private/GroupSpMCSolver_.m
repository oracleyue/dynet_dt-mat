function [w, gtopol, gprob, nvar] = GroupSpMCSolver_(dy, A, indexLists, optimizerOpt)
% Solve "dy = A wi + lambda ||wi||_0^GS" using MCMC methods.
%
% Inputs:
%     indexLists = {[Nblk, Cexp], NblkSize}, 1x2 cell arrays
% Ouputs:
%     w      : estimation of parameters
%     gtopol : group structure (boolean vector)
%     gprob  : prob of group vectors being 1 or 0
%     nvar   : noise variacne

% Copyright [2017] <oracleyue>
% Last modified on 27 Sep 2018


%
% unpack arguments
%
N = indexLists{1}(1);
C = indexLists{1}(2);
blksize = indexLists{2};    % row vector
[Ny, Nw] = size(A);
grp_size = C*blksize;

confid = optimizerOpt.confid;
hypara = optimizerOpt.hypara;
mcpara = optimizerOpt.mcpara;


%
% call solvers to get group sparse structure of w
%

[sG, nvar, gamma] = groupSpMC(dy, A, grp_size, ...
                          'hypara', hypara, 'mcpara', mcpara);
gtopol = ones(N,1);  % indicate group structure (0 or 1)
gprob = sG;   % confid in prob to have sG(i) being 1 or 0
for index = 1:length(sG)
    if sG(index) < .5
        gtopol(index) = 0;
        gprob(index) = 1 - gprob(index);
    end
end


%
% expand group sparsity for w and calculate w
%

grp_midx = zeros(N, 2);
grp_midx = [cumsum([1 grp_size(1:end-1)]); cumsum(grp_size)]';

% prune rows in dy and A that correspond to zero groups
zero_cols = [];
for index = 1:length(gtopol)
    if ~gtopol(index)
        zero_grp = grp_midx(index,:);
        zero_cols = [zero_cols zero_grp(1):zero_grp(2)];
    end
end
A(:, zero_cols) = [];

% calculate truncted w
w_tr = A \ dy;

% lift w_tr to w
w = zeros(Nw,1);

for index = 1:length(gtopol)
    if gtopol(index)
        one_grp = grp_midx(index,:);
        one_size = one_grp(2) - one_grp(1) +1;
        w(one_grp(1):one_grp(2)) = w_tr(1:one_size);
        w_tr(1:one_size) = [];
    end
end

% for debugging
if ~isempty(w_tr)
    error('errors happen in lifting w!')
end
