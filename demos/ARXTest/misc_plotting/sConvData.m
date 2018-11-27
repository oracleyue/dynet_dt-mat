solverNameList = {'GIRL1', 'GSBL', 'GSMC'};
nType = 4;  % length(SNRList) or length(pList);
            % SNRList = [0 10 20 40]; pList = [5 10 15 20];
% type = 'SNR';
type = 'Node';

for indexSolver = 1:length(solverNameList)
    switch type
      case 'SNR'
        switch indexSolver
          case 1 % solver: gil1
            load('dynet_gil1_snr.mat')
          case 2 % solver: gsbl
            load('dynet_gsbl_snr.mat')
          case 3 % solver: gspmc
            load('dynet_gspmc_snr.mat')
        end
      case 'Node'
        switch indexSolver
          case 1 % solver: gil1
            load('dynet_gil1_p.mat')
          case 2 % solver: gsbl
            load('dynet_gsbl_p.mat')
          case 3 % solver: gspmc
            load('dynet_gspmc_p.mat')
        end
    end

    for indexType = 1:nType
        for indexData = 1:size(perf_cells,1)
            mTP(indexData,indexType) = perf_cells{indexData,indexType}(1);
            mFP(indexData,indexType) = perf_cells{indexData,indexType}(2);
            mFN(indexData,indexType) = perf_cells{indexData,indexType}(3);
            mTN(indexData,indexType) = perf_cells{indexData,indexType}(4);
            mPrec(indexData,indexType) = perf_cells{indexData,indexType}(5);
            mTPR(indexData,indexType)  = perf_cells{indexData,indexType}(6);
        end
    end

    cTP{indexSolver} = mTP;
    cFP{indexSolver} = mFP;
    cFN{indexSolver} = mFN;
    cTN{indexSolver} = mTN;
    cPrec{indexSolver} = mPrec*100;  % unit: percentage
    cTPR{indexSolver}  = mTPR *100;
end

switch type
  case 'SNR'
    save('convdata_snr.mat', 'cTP', 'cFP', 'cFN', 'cTN', 'cPrec', 'cTPR');
  case 'Node'
    save('convdata_p.mat', 'cTP', 'cFP', 'cFN', 'cTN', 'cPrec', 'cTPR');
end