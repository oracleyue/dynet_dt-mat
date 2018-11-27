% This script computes the mean and std error of Prec/TPR of
% inference results.

% Copyright [2017] <oracleyue>
% Last modified on 05 Sep 2018


% env setting
solverNameList = {'GIRL1', 'GSBL', 'GSMC'};
resPath = './';

% compute mean and std
% row 1: GIRL1; row 2: GSBL

% task: mSNRs
% info: p = 10; SNRList = [0 10 20 40];
load('convdata_snr.mat');
meanPrec_snr = [mean(cPrec{1}); mean(cPrec{2}); mean(cPrec{3})];
stdPrec_snr  = [std(cPrec{1});  std(cPrec{2});  std(cPrec{3}) ];
meanTPR_snr  = [mean(cTPR{1});  mean(cTPR{2});  mean(cTPR{3}) ];
stdTPR_snr   = [std(cTPR{1});  std(cTPR{2});    std(cTPR{3})  ];

% % compute SEM
% nSize = size(cPrec{1},1);
% semPrec_snr = stdPrec_snr / sqrt(nSize);
% semTPR_snr  = stdTPR_snr / sqrt(nSize);

% task: mNodes
% info: SNR = 10; NodeList = [5 10 15 20];
load('convdata_p.mat');
meanPrec_p = [mean(cPrec{1}); mean(cPrec{2}); nanmean(cPrec{3})];
stdPrec_p  = [std(cPrec{1});  std(cPrec{2});  nanstd(cPrec{3}) ];
meanTPR_p  = [mean(cTPR{1});  mean(cTPR{2});  mean(cTPR{3}) ];
stdTPR_p   = [std(cTPR{1});   std(cTPR{2});   std(cTPR{3})  ];

% % compute SEM
% nSize = size(cPrec{1},1);
% semPrec_p = stdPrec_p / sqrt(nSize);
% semTPR_p  = stdTPR_p / sqrt(nSize);

% display
meanPrec_snr, stdPrec_snr, meanTPR_snr, stdTPR_snr
meanPrec_p, stdPrec_p, meanTPR_p, stdTPR_p

% % display
% meanPrec_snr, semPrec_snr, meanTPR_snr, semTPR_snr
% meanPrec_p, semPrec_p, meanTPR_p, semTPR_p