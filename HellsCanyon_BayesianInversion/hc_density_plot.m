function hc_density_plot(x, y, varargin)
% HC_DENSITY_PLOT  2-D posterior density as a heat map.
%
%   hc_density_plot(x, y)
%   hc_density_plot(x, y, 'nbins', 100, 'scale', 'sqrt')
%
% Replaces scatter plots of MCMC samples.  With ~1e6 post-burn-in samples a
% scatter saturates completely: every pixel inside the support is filled, so
% the plot shows the EXTENT of the posterior but nothing about where the
% probability mass actually sits, and the marker alpha needed to avoid that
% makes the tails invisible.  A 2-D histogram shows both at once.
%
% This also matches Gallen & Fernandez-Blanco's plot_basic_results, which
% presents the joint posteriors as 2-D histograms rather than point clouds.
%
% Options:
%   'nbins'    scalar or [nx ny], default 90
%   'scale'    'linear' (default) | 'sqrt' | 'log'
%              'sqrt'/'log' compress a sharply peaked density so the tails
%              stay visible; 'linear' is the honest default.
%   'cmap'     colormap name or matrix, default 'parula'
%   'overlay'  [xo yo] point drawn on top (e.g. the MAP), default []
%   'label'    true (default) prints the overlay's coordinates beside the
%              marker; false suppresses; or pass a char/string for custom
%              text.  Reading values straight off a colour axis is
%              guesswork, so the numbers go on the plot.
%   'labelfmt' sprintf format for the auto label, default '(%.4g, %.4g)'
%
% See also: plot_hc_results, hc_corner_plot

p = inputParser;
addParameter(p, 'nbins',    90);
addParameter(p, 'scale',    'linear', @ischar);
addParameter(p, 'cmap',     'parula');
addParameter(p, 'overlay',  []);
addParameter(p, 'label',    true);
addParameter(p, 'labelfmt', '(%.4g, %.4g)');
parse(p, varargin{:});
o = p.Results;

x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);

if isempty(x)
    text(0.5, 0.5, 'no finite samples', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center');
    return;
end

nb = o.nbins;
if isscalar(nb), nb = [nb nb]; end

[N, xe, ye] = histcounts2(x, y, nb);

% Bin CENTRES, so the image is not offset by half a bin against the axes.
xc = xe(1:end-1) + diff(xe)/2;
yc = ye(1:end-1) + diff(ye)/2;

switch lower(o.scale)
    case 'sqrt', C = sqrt(N);
    case 'log',  C = log10(N + 1);
    otherwise,   C = N;
end

% Empty bins render as background rather than as "lowest density", so the
% support of the posterior reads clearly.
C(N == 0) = NaN;

imagesc(xc, yc, C', 'AlphaData', ~isnan(C'));
set(gca, 'YDir', 'normal', 'Color', 'w');
colormap(gca, o.cmap);

if ~isempty(o.overlay)
    hold on
    plot(o.overlay(1), o.overlay(2), 'p', 'MarkerSize', 15, ...
         'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

    if ~isequal(o.label, false)
        if ischar(o.label) || isstring(o.label)
            txt = char(o.label);
        else
            txt = sprintf(o.labelfmt, o.overlay(1), o.overlay(2));
        end
        % Offset slightly so the box does not sit on the marker, and give
        % it a solid background so it stays readable over dense bins.
        xr = diff(xlim); yr = diff(ylim);
        text(o.overlay(1) + 0.03*xr, o.overlay(2) + 0.04*yr, txt, ...
             'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', ...
             'BackgroundColor', 'w', 'EdgeColor', [0.4 0.4 0.4], ...
             'Margin', 1, 'Clipping', 'on');
    end
end

axis tight
box on

end
