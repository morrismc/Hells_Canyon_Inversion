function fig = hc_corner_plot(params_post, params_map, param_names, ...
                              param_scale, varargin)
% HC_CORNER_PLOT  Matrix plot of the joint posterior, after Gallen &
% Fernandez-Blanco (2021) Figure 6.
%
%   fig = hc_corner_plot(params_post, params_map, param_names, param_scale)
%   fig = hc_corner_plot(..., 'scale', 'sqrt', 'nbins', 60)
%
% Diagonal   : marginal posterior for each parameter, with the MAP marked
%              and its value annotated (as Gallen prints the MAP value in
%              the corner of each diagonal panel).
% Lower left : bivariate posterior density as a 2-D histogram, showing the
%              covariance between every pair of parameters.
%
% This is the single most informative view of a correlated posterior: the
% marginals alone hide exactly the trade-offs that make a chain mix badly
% (here ksn-n, t_capture-U_post, n-U_grad), and 1e6 samples overplot into a
% solid block in a scatter.
%
% Options:
%   'scale'  'linear' | 'sqrt' (default) | 'log'  -- density colour scaling
%   'nbins'  bins per axis, default 55
%   'maxn'   thin to at most this many samples for speed, default 3e5
%   'cmap'   colormap, default 'parula'
%
% See also: hc_density_plot, plot_hc_results

p = inputParser;
addParameter(p, 'scale', 'sqrt', @ischar);
addParameter(p, 'nbins', 55);
addParameter(p, 'maxn',  3e5);
addParameter(p, 'cmap',  'parula');
% fontscale multiplies EVERY font size, so one number takes the figure
% from screen to poster. >1.3 also thins the ticks, because at poster
% sizes the default tick labels collide across 7x7 panels.
addParameter(p, 'fontscale', 1, @isscalar);
addParameter(p, 'figpos', [60 60 1500 1150]);
parse(p, varargin{:});
o = p.Results;

np = size(params_post, 2);
if nargin < 4 || isempty(param_scale), param_scale = ones(1, np); end
if nargin < 3 || isempty(param_names)
    param_names = arrayfun(@(k) sprintf('p_%d', k), 1:np, 'Uni', false);
end

% Thin for speed.  Evenly spaced rather than random so the thinned set
% still spans the whole chain.
P = params_post;
if size(P,1) > o.maxn
    P = P(round(linspace(1, size(P,1), o.maxn)), :);
end

% Apply display scaling once.
P    = P    .* param_scale(:)';
pmap = params_map(:)' .* param_scale(:)';

fs = o.fontscale;
nt = 3;                      % ticks per axis when scaled up
fig = figure('Position', o.figpos, 'Color', 'w');

for i = 1:np          % row    -> y variable
    for j = 1:i       % column -> x variable
        ax = subplot(np, np, (i-1)*np + j);

        if i == j
            % ---- marginal ----
            histogram(P(:,i), 60, 'Normalization', 'pdf', ...
                      'FaceColor', [0.85 0.82 0.70], ...
                      'EdgeColor', 'none');
            hold on
            yl = ylim;
            plot([pmap(i) pmap(i)], yl, 'r-', 'LineWidth', 1.5*fs);
            ylim(yl);
            % MAP value annotated in the panel corner, as in Gallen Fig 6.
            text(0.96, 0.90, sprintf('%.4g', pmap(i)), ...
                 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'FontSize', 8*fs, 'FontWeight', 'bold', ...
                 'BackgroundColor', 'w', 'Margin', 1);
            set(ax, 'YTick', []);
            title(param_names{i}, 'FontSize', 9*fs, 'Interpreter', 'tex');
        else
            % ---- bivariate density ----
            hc_density_plot(P(:,j), P(:,i), ...
                'nbins', o.nbins, 'scale', o.scale, 'cmap', o.cmap, ...
                'overlay', [pmap(j) pmap(i)], 'label', false, ...
                'markersize', 15*fs);
        end

        set(ax, 'FontSize', 7*fs, 'LineWidth', max(0.5, 0.5*fs));

        % Thin the ticks BEFORE blanking any labels: setting XTick puts the
        % tick labels back into automatic mode, so doing it afterwards would
        % undo the blanking on the interior panels.
        if fs > 1.3
            % Fewer, rounder ticks so poster-size labels do not overlap.
            xl = xlim(ax);
            set(ax, 'XTick', linspace(xl(1), xl(2), nt));
            xtickformat(ax, '%.3g');
            if i ~= j
                ylv = ylim(ax);
                set(ax, 'YTick', linspace(ylv(1), ylv(2), nt));
                ytickformat(ax, '%.3g');
            end
        end

        % Label only the outer edges, or the grid becomes unreadable.
        if i == np
            xlabel(param_names{j}, 'FontSize', 8*fs, 'Interpreter', 'tex');
        else
            set(ax, 'XTickLabel', []);
        end
        if j == 1 && i > 1
            ylabel(param_names{i}, 'FontSize', 8*fs, 'Interpreter', 'tex');
        elseif j > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

sgtitle({'Posterior probability matrix', ...
         'diagonal: marginals (MAP in red);  lower: bivariate density'}, ...
        'FontSize', 13*fs, 'Interpreter', 'none');

end
