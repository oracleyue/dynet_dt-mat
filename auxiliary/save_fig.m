function save_fig(fig_name, fig_hl, fig_pos, units)
% SAVE_FIG Save figures in PDF, using the dimension of figure shown,
%
% Copyright [2017] <oracleyue>


if nargin == 1
    fig_hl = gcf;
elseif nargin == 2
    ;
elseif nargin == 3
    units = 'pixels';
    warning('The default unit is pixels. Use the fourth argument to change units.')
    if length(fig_pos) ~= 4
        error('figure.Position should be 1x4 vector in double.')
    end
    set(fig_hl,'Units', units, 'Position', fig_pos);
else
    if length(fig_pos) ~= 4
        error('figure.Position should be 1x4 vector in double.')
    end
    set(fig_hl,'Units', units, 'Position', fig_pos);
end

match_status = regexp(fig_name, '\.pdf$', 'once');
if isempty(match_status)
    fig_name = [fig_name, '.pdf'];
end


set(fig_hl,'Units','Inches');
pos = get(fig_hl,'Position');
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, fig_name, '-dpdf', '-r0')
