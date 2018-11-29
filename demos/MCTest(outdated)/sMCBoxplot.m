% Bar plots of Type-I/II errors using different methods in the Monte
% Carlo study.

% Copyright [2017] <oracleyue>



% clear all; close all


% task setting
taskSelected = 'multiple SNRs';
% taskSelected = 'multiple Nodes';
solverNameList = {'GIRL1', 'GSBL'};

% read file and organize into "mTP, mFP, mFN" (m-: matrix; boxchild_hl-: cell)
% info: row - datasets; col - different SNR or nodes
switch taskSelected
  case 'multiple SNRs'
    resPath = '../../Results/resSNRs/';

    if ~exist('p') || ~exist('SNRList')
        p = 10;                     % dim of y
        SNRList = [80 60 40 20 0];  % signal to noise ratio
    end

    for indexSolver = 1:length(solverNameList)
        mTP = []; mFP = []; mFN = []; mTN = [];
        for index = 1:length(SNRList)
            SNR = SNRList(index);
            txtName = [resPath 'dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_' ...
                       solverNameList{indexSolver} '.txt'];
            results = dlmread(txtName, '\t', 1, 0);
            mTP(:,index) = results(:,1);
            mFP(:,index) = results(:,2);
            mFN(:,index) = results(:,3);
            mTN(:,index) = results(:,4);
        end
        cTP{indexSolver} = mTP;
        cFP{indexSolver} = mFP;
        cFN{indexSolver} = mFN;
        cTN{indexSolver} = mTN;
    end


  case 'multiple Nodes'
    resPath = '../../Results/resNodes/';

    if ~exist('SNR') || ~exist('NodeList')
        SNR = 40;
        NodeList = [10 20 40 80];
    end

    for indexSolver = 1:length(solverNameList)
        mTP = []; mFP = []; mFN = []; mTN = [];
        for index = 1:length(SNRList)
            p = NodeList(index);
            txtName = [resPath 'dynet_p' num2str(p) '_SNR' num2str(SNR) 'dB_' ...
                       solverNameList{indexSolver} '.txt'];
            results = dlmread(txtName, '\t', 1, 0);
            mTP(:,index) = results(:,1);
            mFP(:,index) = results(:,2);
            mFN(:,index) = results(:,3);
            mTN(:,index) = results(:,4);
        end
        cTP{indexSolver} = mTP;
        cFP{indexSolver} = mFP;
        cFN{indexSolver} = mFN;
        cTN{indexSolver} = mTN;
    end
end
for indexSolver = 1:length(solverNameList)
    mTP = cTP{indexSolver};
    mFP = cFP{indexSolver};
    mFN = cFN{indexSolver};

    mPrec = mTP ./ (mTP + mFP);
    mTPR  = mTP ./ (mTP + mFN);
    typeIError  = (ones(size(mPrec)) - mPrec)*100; % unit: percentage
    typeIIError = (ones(size(mTPR)) - mTPR)*100;
    cTypeIError{indexSolver}  = typeIError;
    cTypeIIError{indexSolver} = typeIIError;
end


% Multiple Grouped Boxplots

% xlabel for figures
switch taskSelected
  case 'multiple SNR'
    xAxisNames = strread(num2str(SNRList), '%s')';
    xLabelName = 'SNR (dB)';
  case 'multiple Nodes'
    xAxisNames = strread(num2str(NodeList), '%s')';
    xLabelName = '#Nodes';
end

% group data
dataTypeIError = []; dataTypeIIError = [];
endIndex = length(xAxisNames)*length(solverNameList);
for indexSolver = 1:length(solverNameList)
    dataTypeIError= [dataTypeIError cTypeIError{indexSolver}];
    dataTypeIIError= [dataTypeIIError cTypeIIError{indexSolver}];
end
indexReshape = [1:1:endIndex/2; endIndex/2+1:1:endIndex];  % specific for length(solverNameList) == 2
indexReshape = reshape(indexReshape, 1, []);
dataTypeIError = dataTypeIError(:,indexReshape);
dataTypeIIError = dataTypeIIError(:,indexReshape);

% option variable for boxplot
innerSpace = 1;
pairSpace = innerSpace + 2.5;
positions = (1:pairSpace:1+pairSpace*(length(xAxisNames)-1))';
positions = [positions positions+innerSpace];
positions = reshape(positions', 1, []);
xtickPositions = mean(reshape(positions, 2, []));

% figure name
if ~exist('figNamePrefix')
    figName = [resPath 'out_err'];
end

% setup colormap
colorTable = colormap('lines');
colorSelected = colorTable(1:2, :);
color = repmat(colorSelected, length(xAxisNames), 1);
close

% Figure 1: type-I error
fig_hl = figure('visible', 'off');
set(fig_hl, 'units', 'inches', 'position',[9.4306 11.2778 7.4722 2.2917])

bp_hl = boxplot(dataTypeIError, 'position', positions);

xlabel(xLabelName);
ylabel('Type-I Error');
set(gca,'xtick', xtickPositions);
set(gca,'xticklabel', xAxisNames);

% color = repmat('br', 1,4);
box_hl = findobj(gca,'Tag','Box');
for j=1:length(box_hl)
   patch(get(box_hl(j),'XData'),get(box_hl(j),'YData'), color(j,:),'FaceAlpha',.8);
   set(bp_hl(6,j), 'color', [0 0 0]);  %'linewidth'
end

boxch_hl = get(gca, 'Children');
lg_hl = legend(boxch_hl(1:2), 'GIRL1', 'GSBL' );
% set(lg_hl, 'location', 'northwest');

save_fig([figName '1.pdf'], fig_hl);


% Figure 2: type-II error
fig_hl = figure('visible', 'off');
set(fig_hl, 'units', 'inches', 'position',[9.4306 11.2778 7.4722 2.2917])

bp_hl = boxplot(dataTypeIIError, 'position', positions);

xlabel(xLabelName);
ylabel('Type-II Error');
set(gca,'xtick', xtickPositions);
set(gca,'xticklabel', xAxisNames);

% color = repmat('br', 1,4);
box_hl = findobj(gca,'Tag','Box');
for j=1:length(box_hl)
   patch(get(box_hl(j),'XData'),get(box_hl(j),'YData'),color(j,:),'FaceAlpha',.8);
   set(bp_hl(6,j), 'color', [0 0 0]);
end

boxch_hl = get(gca, 'Children');
lg_hl = legend(boxch_hl(1:2), 'GIRL1', 'GSBL' );
% set(lg_hl, 'location', 'northwest');

save_fig([figName '2.pdf'], fig_hl);



% ----------------------------------------------------------------
% % ARCHIVE: Single Grouped Boxplots
% ----------------------------------------------------------------

% % figure 1: type-I error
% fig_hl = figure;
% bh = boxplot(typeIError, xAxisNames);
% for i=1:size(bh,1) % <- # graphics handles/x
%     for j = 1:size(bh,2)
%         set(bh(i,j),'linewidth',2);
%     end
% end
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 1;
% ax.Color = 'none'; % OR set(gca, 'Color', 'none'); % transparent bg color
% xlabel('SNR (dB)', 'FontSize', 16)
% ylabel('Type-I Error (%)', 'FontSize', 16)

% set(fig_hl,'Units','Inches', 'Position', [7.0208 5.7396 5.9375 4.3958]);
% export_fig -transparent snr_out1.pdf


% % figure 2: type-II error
% fig_hl = figure;
% bh = boxplot(typeIIError, xAxisNames);
% for i=1:size(bh,1) % <- # graphics handles/x
%     for j = 1:size(bh,2)
%         set(bh(i,j),'linewidth',2);
%     end
% end
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 1;
% ax.Color = 'none'; % OR set(gca, 'Color', 'none'); % transparent bg color
% xlabel('SNR (dB)', 'FontSize', 16)
% ylabel('Type-II Error (%)', 'FontSize', 16)

% set(fig_hl,'Units','Inches', 'Position', [7.0208 5.7396 5.9375 4.3958]);
% export_fig -transparent snr_out2.pdf
