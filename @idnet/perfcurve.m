function fig_hl = perfcurve(net, boolNets, curve_type, varParaVector)
% PERFCURVE is to draw the performance curves, e.g. P-R, ROC, etc.
%
% Input:
%   boolNets: the cell of bool networks returned by `multi-run()`
%   curve_type: character array, can be:
%       - 'prec'  : Precision-Recall curve
%       - 'roc'   : ROC curve
%       - 'err'   : Type-I and Type-II error plots
%       - 'all'   : the above three curves
%
% Output:
%   fig_hl: the handle of the figure
%
%
% Copyright [2017] <oracleyue>
%


switch curve_type
  case 'prec'
    fig_hl = prec_recall(boolNets, net.GroundTruth);

  case 'roc'
    fig_hl = roc(boolNets, net.GroundTruth);

  case 'err'
    fig_hl = topolo_err(boolNets, net.GroundTruth, varParaVector);

  case 'all'
    fig_hl1 = prec_recall(boolNets, net.GroundTruth);
    fig_hl2 = roc(boolNets, net.GroundTruth);
    fig_hl3 = topolo_err(boolNets, net.GroundTruth, varParaVector);
    fig_hl = {fig_hl1, fig_hl2, fig_hl3};

  otherwise
    error('curve_type is invalid! e.g. prec, roc, err.')

end
