function fig_hl = prec_recall(boolNets, groundTruth)
% PREC_RECALL is to plot the precision-recall curve in terms of the
% given variable (e.g. model orders, labmdas).
%   P-R curve: TP/(TP+FN) (Recall) vs. TP/(TP+FP) (Precision)
%   ROC curve: FP/(FP+TN) (FPR) vs. TP/(TP+FN) (TPR)
%
% Input:
%   boolNets: 1xN cell
%   groundTruth: adjacency matrix
%
% Copyright [2017] <oracleyue>


% Extract dimensions
p = size(groundTruth,1);
N = length(boolNets);

% Statistics of inferred results
boolNetSum = zeros(p,p);
TP = zeros(N,1); FP = zeros(N,1); FN = zeros(N,1);
for i = 1:N
    boolNetSum = boolNetSum + boolNets{i}.Q';
    % TP(i) = sum(sum( groundTruth &  boolNets{i}.Q'));
    % FP(i) = sum(sum(~groundTruth &  boolNets{i}.Q'));
    % FN(i) = sum(sum( groundTruth & ~boolNets{i}.Q'));
end
% TypeIError  = 1 - TP./(TP + FP);
% TypeIIError = 1 - TP./(TP + FN);

% Calculate variables for Precision-Recall Curves
boundaryPoints = unique(boolNetSum);
Npr = length(boundaryPoints);
TPpc = zeros(Npr,1); TNpc = zeros(Npr,1);
FPpc = zeros(Npr,1); FNpc = zeros(Npr,1);
for i = 1:Npr
    ptVal = boundaryPoints(i);
    newBoolNetSum = (boolNetSum >= ptVal);
    TPpc(i) = sum(sum( groundTruth &  newBoolNetSum));
    TNpc(i) = sum(sum(~groundTruth & ~newBoolNetSum));
    FPpc(i) = sum(sum(~groundTruth &  newBoolNetSum));
    FNpc(i) = sum(sum( groundTruth & ~newBoolNetSum));
end
% For Precision-Recall curve
Recall = TPpc./(TPpc + FNpc);
Precision = TPpc./(TPpc + FPpc);
pts = [Recall, Precision];
% PRpts = StairReshape(pts, 'pr');
ptRaw = [0 1; pts];
ptRaw = flipud(sortrows(ptRaw));
[~,indexSort] = sort(ptRaw(:,1));
PRpts = ptRaw(indexSort,:);

% % For ROC curve
% FPR = FPpc./(FPpc+TNpc);
% TPR = TPpc./(TPpc+FNpc);
% pts = [FPR, TPR];
% ROCpts = StairReshape(pts, 'roc');


% Figures Plotting

% % Fig 1. ModelOrders vs. Type-I/II Error
% %   Type I  Error = 1 - TP/(TP + NP)
% %   Type II Error = 1 - TP/(TP + FN)
% figure('Visible', 'on');
% plot(TypeIError, '--o');
% hold on;
% plot(TypeIIError, '--s');
% hold off;
% legend('Type I Error: NP/(TP+NP)', 'Type II Error: FN/(TP+FN)');
% xlabel('Orders of Parametric Models');

% Fig 2. Precision-Recall Curve in terms of model orders
fig_hl = figure('Visible', 'on');
plot(PRpts(:,1), PRpts(:,2), '-s');
hold on;
plot([min(PRpts(:,1)) max(PRpts(:,1))], ...
     [-min(PRpts(:,1))+1 -max(PRpts(:,1))+1], '-.', ...
     'color', [.5 .5 .5]);
hold off;
xlim([0 1]); ylim([0 1]);
title('Precision-Recall Curve for Model Orders');
xlabel('Recall = TP/(TP+FN)');
ylabel('Precision = TP/(TP+FP)');
legend('P-R curve', 'y=1-x');

% % Fig 3. Precision-Recall Curve in terms of model orders
% figure('Visible', 'on');
% plot(ROCpts(:,1), ROCpts(:,2), '-bo');
% hold on;
% plot([min(ROCpts(:,1)) max(ROCpts(:,1))], [min(ROCpts(:,1)) max(ROCpts(:,1))], '-.');
% hold off;
% title('ROC Curve for Model Orders');
% xlabel('FPR = FP/(FP+TN)');
% ylabel('TPR = TP/(TP+FN)');

end




% % ------- Local Functions -------
% function ptVec = StairReshape(ptRaw, type)
% % sorting data points in PR/ROC curves and insert more to
% % show stair-like shapes
%
% switch type
%   case 'roc'
%     % sort in ascending order
%     ptRaw  = sortrows(ptRaw);
%     % additional points to have stair shapes
%     ptAdd = [ptRaw(2:end, 1) ptRaw(1:end-1, 2)];
%     % merge data points
%     ptVec = [ptRaw; ptAdd];
%     ptVec = sortrows(ptVec);
%
%   case 'pr'
%     ptRaw = flipud(sortrows(ptRaw));
%     [~,indexSort] = sort(ptRaw(:,1));
%     ptRaw = ptRaw(indexSort,:);
%
%     ptAdd = [ptRaw(2:end, 1) ptRaw(1:end-1, 2)];
%
%     ptVec = [];
%     for i = 1:size(ptRaw,1)-1
%         ptVec = [ptVec; ptRaw(i,:); ptAdd(i,:)];
%     end
% end
%
% end
