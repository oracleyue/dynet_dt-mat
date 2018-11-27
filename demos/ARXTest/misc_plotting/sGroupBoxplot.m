% Grouped boxplots of Prec/TPR of results inferred by GIRL1/GSBL
% for ARX dynet models.

% Copyright [2018] <oracleyue>
% Last modified on 04 Sep 2018


clear all; close all

% task setting
% taskSelected = 'mSNRs';
taskSelected = 'mNodes';
solverNameList = {'GIRL1', 'GSBL', 'GSMC'};

resPath = './';
switch taskSelected
  case 'mSNRs'
    p = 10;                     % dim of y
    SNRList = [0 10 20 40];     % signal to noise ratio
    load('convdata_snr.mat');
  case 'mNodes'
    SNR = 10;
    NodeList = [5 10 15 20];
    load('convdata_p.mat');
end


% Multiple Grouped Boxplots

% xlabel for figures
switch taskSelected
  case 'mSNRs'
    xAxisNames = strread(num2str(SNRList), '%s')';
    xLabelName = 'SNR (dB)';
  case 'mNodes'
    xAxisNames = strread(num2str(NodeList), '%s')';
    xLabelName = '#nodes';
end

% group data
dataPrec = []; dataTPR = [];
endIndex = length(xAxisNames)*length(solverNameList);
for indexSolver = 1:length(solverNameList)
    dataPrec= [dataPrec cPrec{indexSolver}];
    dataTPR= [dataTPR cTPR{indexSolver}];
end
% % specific for length(solverNameList) == 2
% indexReshape = [1:1:endIndex/2; endIndex/2+1:1:endIndex];
% specific for length(solverNameList) == 3
indexReshape = [1:1:endIndex/3;
                endIndex/3+1:1:endIndex/3*2;
                endIndex/3*2+1:1:endIndex];
indexReshape = reshape(indexReshape, 1, []);
dataPrec = dataPrec(:,indexReshape);
dataTPR = dataTPR(:,indexReshape);

% option variable for boxplot
innerSpace = .5;
% % specific for length(solverNameList) == 2
% pairSpace = innerSpace + 1;
% positions = (1:pairSpace:1+pairSpace*(length(xAxisNames)-1))';
% positions = [positions positions+innerSpace];
% specific for length(solverNameList) == 3
pairSpace = innerSpace*2 + 1;
positions = (1:pairSpace:1+pairSpace*(length(xAxisNames)-1))';
positions = [positions positions+innerSpace positions+2*innerSpace];
xTickPositions = reshape(positions', 1, []);
groupPositions = mean(positions, 2);

% figure name
if ~exist('figNamePrefix')
    figName = [resPath taskSelected];
end

% setup colormap
colorTable = colormap('lines');
colorSelected = colorTable(1:3, :);
color = repmat(colorSelected, length(xAxisNames), 1);
close

% Figure 1: type-I error
% fig_hl = figure('visible', 'off');
fig_hl = figure(1);
set(fig_hl, 'units', 'inches', 'position',[9.4306 11.2778 7.4722 2.2917])

bp_hl = boxplot(dataPrec, 'position', xTickPositions);

xlabel(xLabelName);
ylabel('Prec (%)');
set(gca,'xtick', groupPositions);
set(gca,'xticklabel', xAxisNames);

% color = repmat('br', 1,4);
box_hl = findobj(gca,'Tag','Box');
for j=1:length(box_hl)
   patch(get(box_hl(j),'XData'),get(box_hl(j),'YData'), color(j,:),...
         'FaceAlpha',.8);
   set(bp_hl(6,j), 'color', [0 0 0]);  %'linewidth'
end

boxch_hl = get(gca, 'Children');
lg_hl = legend(boxch_hl(1:3), solverNameList);
set(lg_hl, 'location', 'southeast');

save_fig([figName '_Prec.pdf'], fig_hl);


% Figure 2: type-II error
% fig_hl = figure('visible', 'off');
fig_hl = figure(2);
set(fig_hl, 'units', 'inches', 'position',[9.4306 11.2778 7.4722 2.2917])

bp_hl = boxplot(dataTPR, 'position', xTickPositions);

xlabel(xLabelName);
ylabel('TPR (%)');
set(gca,'xtick', groupPositions);
set(gca,'xticklabel', xAxisNames);

% color = repmat('br', 1,4);
box_hl = findobj(gca,'Tag','Box');
for j=1:length(box_hl)
   patch(get(box_hl(j),'XData'),get(box_hl(j),'YData'),color(j,:),...
         'FaceAlpha',.8);
   set(bp_hl(6,j), 'color', [0 0 0]);
end

boxch_hl = get(gca, 'Children');
lg_hl = legend(boxch_hl(1:3), solverNameList);
set(lg_hl, 'location', 'southeast');

save_fig([figName '_TPR.pdf'], fig_hl);
