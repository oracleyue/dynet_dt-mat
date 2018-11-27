function [fit, y, x0] = compare(dsf, ndata, plotFlag)
% To compare model outpu and measured output. Modified from /compare/ in
% sys.id toolbox for dynamic network reconstruction, using the same declaration.
%
% Input:
%   - ndata:    iddata
%               input-output data for system id
%   - plotFlag: bool or string
%               1 or "plot" - plot and export in PDF;
%               0 (default) - not plot
%
% Output:
%   - fit: px1 vector
%          NRMSE fitness value, in percentage, using Inf prediction horzion.
%   - y:   px1 cell
%          Model response. Each element is the model response of an MISO sys.
%   - x0:  px1 cell
%          Initial conditions used to compute system response. Each element
%          corresponds to an MISO sys.
% Usages:
%   - compare(..)
%     Plot comparation results and save them in PDF under './Results/compare_results/'.
%   - fit = compare(..)
%   - [fit, y] = compare(..)
%   - [fit, y, x0] = compare(..)
%
%
% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 05 Dec 2016


% Input argumenets parsing
numArgumentIn = nargin;
numArgumentOut = nargout;

if ~numArgumentOut
    plotFlag = 1;
end

if numArgumentIn < 2
    error('Not enough input arguments!')
elseif numArgumentIn == 2 && numArgumentOut
    plotFlag = 0;
    disp 'Im here'
elseif numArgumentIn == 3
    if strcmp(plotFlag, 'plot') || strcmp(plotFlag, 'Plot')
        plotFlag = 1;
    else
        plotFlag = 0;
    end
elseif numArgumentIn > 3
    error('Too many input arguments!')
end


% Distill model information
p = size(dsf.Q, 1);
m = size(dsf.P, 2);
Ts = dsf.Ts;

% Path to save "compare" plots
if ~exist('Results', 'dir')
    mkdir('Results', 'dir')
end
if ~exist('Results/compare_results', 'dir')
    mkdir('Results/compare_results');
end
respath = './Results/compare_results/';

% Variables for collecting values in parfor
fit = zeros(p,1);  % name follows /compare/ in sys.id. toolbox

% Comparing process
type = 'Q';
switch type
  case 'Q'
    parfor index = 1:p
        sys = dsf.Q(index,:); sys(index) = [];

        output = ndata.y(:,index);
        input = ndata.y; input(:,index) = [];
        data = iddata(output, input, Ts);

        [modelResponse, fitScore, InitCondition] = compare(data, sys);
        fit(index) = fitScore;
        if numArgumentOut == 2
            y{index} = modelResponse;
        elseif numArgumentOut == 3
            y{index} = modelResponse;
            x0{index} = InitCondition;
        end
    end

  case 'QP'
    error('The function is NOT yet available!');

end

% Set negative values (e.g. -Inf) to be zeros
% fit(find(fit < 0)) = 0;

% Plotting the model comaprison
if plotFlag

    % Default figures settings
    set(0, 'DefaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 2);

    % Plotting all channels
    ndimY = size(ndata.y, 2);
    for k = 1:ndimY
        dataYk = ndata.y(:,k);
        dataSimYk = y{k}.y;
        fitMk = fit(k);

        Ts = ndata.Ts;
        Nsamples = length(dataYk);
        time = 0:Ts:Ts*(Nsamples-1);

        fig_hl = figure('visible', 'off');
        subplot(2,1,1);
        plot(time, dataYk, time, dataSimYk);
        title('Model Comparison')
        xlabel(['time (Ts = ' num2str(Ts) ')']);
        ylabel(['y' num2str(k)]);
        legend('data', ['model (' 'fit: ' num2str(fitMk) '%)']);
        subplot(2,1,2);
        set(gca, 'ColorOrderIndex', 4); hold on;
        plot(time, dataSimYk - dataYk, 'LineWidth', 1); box on;
        xlabel(['time (Ts = ' num2str(Ts) ')']);
        ylabel('prediction errors'); hold off

        % Export the figure in PDF
        figname = [respath 'fitness_y' num2str(k) '.pdf'];
        set(fig_hl,'Units','Inches');
        pos = get(fig_hl,'Position');
        set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(fig_hl, figname, '-dpdf', '-r0')
    end

end
